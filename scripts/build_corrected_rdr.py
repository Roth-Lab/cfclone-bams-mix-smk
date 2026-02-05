from statsmodels.nonparametric.smoothers_lowess import lowess

import numpy as np
import pandas as pd
import pybedtools
import scipy
import statsmodels.formula.api as smf


def main(args):
    df_anno = get_annotations(
        args.black_list_file,
        args.centromere_file,
        args.gc_file,
        args.map_file,
        args.normal_reads_file,
    )

    df = pd.read_csv(args.in_file, converters={"chrom": str}, sep="\t")

    df = pd.merge(df, df_anno, on=["chrom", "start", "end"])

    df = gc_correction(df, min_mappability=args.min_mappability)

    df.to_csv(args.out_file, compression="gzip", index=False, sep="\t")


def get_annotations(
        black_list_file,
        centromere_file,
        gc_file,
        map_file,
        normal_file,
):
    df = pd.read_csv(normal_file, sep="\t")

    df = df[["chrom", "start", "end", "reads"]]

    df = df.rename(columns={"reads": "normal"})

    df["bin"] = df.apply(lambda row: "{chrom}:{start}:{end}".format(**row), axis=1)

    # Blacklist info
    bins_bed = pybedtools.BedTool(list(df[["chrom", "start", "end", "bin"]].values))

    df = df.set_index("bin")

    black_list_bed = pybedtools.BedTool(black_list_file)

    df["black_list"] = False

    for x in bins_bed.intersect(black_list_bed):
        df.loc[x.name, "black_list"] = True

    # Centromere info
    centromere_df = pd.read_csv(centromere_file, sep="\t")

    centromere_bed = pybedtools.BedTool(list(centromere_df.values))

    df["centromere"] = False

    for x in bins_bed.intersect(centromere_bed):
        df.loc[x.name, "centromere"] = True

    df = df.reset_index()

    df = append_bin_info(["gc"], df, gc_file)

    df = append_bin_info(["mappability"], df, map_file)

    df = df.drop(columns="bin")

    return df


def append_bin_info(cols, df, file_name):
    bin_df = pd.read_csv(file_name, sep="\t")

    if "chr" in bin_df.columns:
        bin_df = bin_df.rename(columns={"chr": "chrom"})

    bin_df["bin"] = bin_df.apply(
        lambda row: "{chrom}:{start}:{end}".format(**row.to_dict()), axis=1
    )

    cols = ["bin"] + cols

    return pd.merge(df, bin_df[cols], on="bin", how="left")


def gc_correction(df, min_mappability=0.9):
    scale = df["normal"].sum() / df["reads"].sum()

    ratio = df["reads"] / df["normal"]

    df["rdr"] = ratio * scale

    df["valid"] = (
            ~df["black_list"]
            & ~df["centromere"]
            & (df["gc"] > 0)
            & (df["mappability"] > min_mappability)
    )

    # filtering and sorting
    df_regression = df[df["valid"]].copy()

    df_regression.sort_values(by="gc", inplace=True)

    df_regression = modal_quantile_regression(df_regression, lowess_frac=0.2)

    df["gc_correction"] = np.nan

    df["rdr_cor"] = np.nan

    df.loc[df_regression.index, "gc_correction"] = df_regression["modal_curve"]

    df["rdr_cor"] = df["rdr"] / df["gc_correction"]

    return df


def modal_quantile_regression(df_regression, lowess_frac=0.2, degree=2, knots=(0.38,)):
    """
    Fits a B-spline polynomial curve through the "modal" quantile of the data:
    * Runs quantile regression to fit a B-spline curve for each percentile 10-90
    * Estimates the modal quantile as the quantile where difference in AUC is minimized
    * Uses the curve fit to this modal quantile for normalization

    Parameters:
        df_regression: pandas.DataFrame with at least columns [chr, start, end, rdr, gc]
        lowess_frac: float, fraction of data used to estimate each y-value in Lowess smoothing of AUC curve
        degree: int, degree of polynomial to fit to each section of the B-spline curve
        knots: list of floats, GC values where B-spline polynomial is allowed to change

    Returns:
        pandas.DataFrame with additional columns
            modal_curve: modal curve's predicted # rdr for GC value in this row
            modal_quantile: quantile selected as the mode (should be the same for all bins)
            modal_corrected: corrected read count (i.e., rdr / modal_curve)
    """

    q_range = range(10, 91, 1)
    quantiles = np.array(q_range) / 100
    quantile_names = [str(x) for x in q_range]

    # need at least 3 values to compute the quantiles
    if len(df_regression) < 10 or sum(df_regression["rdr"]) < 100:
        df_regression["modal_quantile"] = None
        df_regression["modal_curve"] = None
        df_regression["modal_corrected"] = None
        return df_regression

    poly_quantile_model = smf.quantreg(
        f"rdr ~ bs(gc, degree={degree}, knots={knots}, include_intercept = True)",
        data=df_regression,
    )

    poly_quantile_fit = [
        poly_quantile_model.fit(q=q, max_iter=10000) for q in quantiles
    ]

    poly_quantile_predict = [
        poly_quantile_fit[i].predict(df_regression) for i in range(len(quantiles))
    ]

    poly_quantile_params = pd.DataFrame()

    for i in range(len(quantiles)):
        df_regression[quantile_names[i]] = poly_quantile_predict[i]
        poly_quantile_params[quantile_names[i]] = poly_quantile_fit[i].params

    # integration and mode selection
    gc_min = df_regression["gc"].quantile(q=0.10)
    gc_max = df_regression["gc"].quantile(q=0.90)

    true_min = df_regression["gc"].min()
    true_max = df_regression["gc"].max()

    poly_quantile_integration = np.zeros(len(quantiles) + 1)

    # form (k+1)-regular knot vector
    repeats = degree + 1

    my_t = np.r_[[true_min] * repeats, knots, [true_max] * repeats]

    for i in range(len(quantiles)):
        # compose params into piecewise polynomial
        params = poly_quantile_params[quantile_names[i]].to_numpy()
        pp = scipy.interpolate.PPoly.from_spline((my_t, params[1:] + params[0], degree))

        # compute integral
        poly_quantile_integration[i + 1] = pp.integrate(gc_min, gc_max)

    # find the modal quantile
    distances = poly_quantile_integration[1:] - poly_quantile_integration[:-1]

    df_dist = pd.DataFrame(
        {
            "quantiles": quantiles,
            "quantile_names": quantile_names,
            "distances": distances,
        }
    )
    dist_max = df_dist["distances"].quantile(q=0.95)
    df_dist_filter = df_dist[df_dist["distances"] < dist_max].copy()
    df_dist_filter["lowess"] = lowess(
        df_dist_filter["distances"],
        df_dist_filter["quantiles"],
        frac=lowess_frac,
        return_sorted=False,
    )

    modal_quantile = df_dist_filter.set_index("quantile_names")["lowess"].idxmin()

    # add values to table
    df_regression["modal_quantile"] = modal_quantile

    df_regression["modal_curve"] = df_regression[modal_quantile]

    return df_regression


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-b", "--black-list-file", required=True)

    parser.add_argument("-c", "--centromere-file", required=True)

    parser.add_argument("-g", "--gc-file", required=True)

    parser.add_argument("-m", "--map-file", required=True)

    parser.add_argument("-n", "--normal-reads-file", required=True)

    parser.add_argument("--min_mappability", default=0.9, type=float)

    cli_args = parser.parse_args()

    main(cli_args)

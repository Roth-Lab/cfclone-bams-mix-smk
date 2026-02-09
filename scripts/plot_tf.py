import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt


plot_settings = {
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
}


def plot(
    df_plot: pd.DataFrame,
    df_patient_tc: pd.DataFrame | None = None,
    figsize: tuple[int, int] | None = None,
    output_path: str | None = None,
    show: bool = False
) -> None:
    with plt.rc_context(plot_settings):
        _plot(
            df_plot=df_plot,
            df_patient_tc=df_patient_tc,
            figsize=figsize,
            output_path=output_path,
            show=show
        )
    
def _plot(
    df_plot: pd.DataFrame,
    df_patient_tc: pd.DataFrame | None = None,
    figsize: tuple[int, int] | None = None,
    output_path: str | None = None,
    show: bool = False
) -> None:
    
    patients = df_plot['patient'].unique()
    
    num_patients = len(patients)
    
    assert num_patients == 1
    
    samples = df_plot['sample'].unique()
    
    num_sample_ids = len(samples)
    
    assert num_sample_ids == 1
    
    coverage_ids = df_plot['coverage'].unique()
    
    num_coverage_ids = len(coverage_ids)
    
    assert num_coverage_ids == 1
    
    num_rows = 1
    
    num_cols = 1
    
    figsize = (num_cols * 6, num_rows * 4)
        
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    
    gs = fig.add_gridspec(num_rows, num_cols)
    
    ax = fig.add_subplot(gs[0, 0])
    
    df_tmp = df_plot.sort_values(by=['proportion'])
    
    ax.scatter(
        x=df_tmp['proportion'], 
        y=df_tmp["mean"]
    )
    
    ax.errorbar(
        x=df_tmp['proportion'],
        y=df_tmp["mean"],
        yerr=[
            df_tmp["mean"] - df_tmp["lower_hdi"],
            df_tmp["upper_hdi"] - df_tmp["mean"],
        ],
        fmt="o",
        ecolor="gray",
        alpha=0.5,
    )
    
    # if df_patient_tc is not None:
    #     df = df_patient_tc.loc[df_patient_tc['sample_id'].isin(['VOA8553P', sample_id])]
    #     ax.errorbar(
    #         x=[0., 1.],
    #         y=df['mean'],
    #         yerr=[
    #             df['mean'] - df['lower_hdi'],
    #             df['upper_hdi'] - df['mean'],
    #         ],
    #         fmt="o",
    #         c='orange',
    #     )
    
    
    sample = df_tmp['sample'].values[0]
    
    coverage = df_tmp['coverage'].values[0]
    
    ax.set_ylabel("Tumour Content estimate")
    
    ax.set_xlabel("Proportion of {}".format(sample))
    
    ax.set_title("Coverage {}X".format(round(coverage, 2)))

    if output_path:
        
        plt.savefig(output_path)

    if show:
        
        plt.show()

    plt.close()
    
    
def main(args):
    
    df = pd.read_csv(args.in_file, sep='\t')
    
    plot(
        df_plot=df,
        show=True,
        output_path=args.out_file,
    )
    

if __name__ == "__main__":
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    
    parser.add_argument('-i', '--in-file', type=str, default="summary.tsv") 
    
    parser.add_argument('-o', '--out-file', type=str, default="lol.png")

    cli_args = parser.parse_args()
    
    main(cli_args)
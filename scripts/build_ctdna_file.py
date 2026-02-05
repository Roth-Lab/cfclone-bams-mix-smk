import pandas as pd


def main(args):
    df = pd.read_csv(args.in_file, sep="\t")

    # df = df.drop(columns=["rdr", "d"])
    df = df.drop(columns=["rdr"])

    # df = df.rename(columns={"rdr_cor": "rdr"})
    df = df.rename(columns={"allele_0_count": "a", "allele_1_count": "b", "rdr_cor": "rdr"})
    
    df = df[df["valid"]]

    df = df[["chrom", "start", "end", "a", "b", "rdr"]].drop_duplicates()
    
    df = df.astype({'a': int, 'b': int})

    df.to_csv(args.out_file, compression="gzip", index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
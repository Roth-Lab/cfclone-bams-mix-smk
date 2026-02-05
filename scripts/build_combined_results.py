import pandas as pd


def main(args):
    baf = pd.read_csv(args.baf_file, sep="\t")

    rdr = pd.read_csv(args.rdr_file, sep="\t")

    df = pd.merge(baf, rdr, on=["chrom", "start", "end"])

    df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--baf-file", required=True)

    parser.add_argument("-r", "--rdr-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)

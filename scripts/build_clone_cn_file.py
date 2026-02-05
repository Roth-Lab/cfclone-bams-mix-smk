import pandas as pd


def main(args):
    clone_df = pd.read_csv(args.clone_file, sep="\t")

    print(clone_df)

    clones = clone_df.loc[clone_df["keep"], "clone_id"].astype(int)

    df = pd.read_csv(args.in_file, sep="\t")

    df = df.rename(columns={"cn_A": "cn_a", "cn_B": "cn_b", "cluster_id": "clone", "beg": "start"})

    df = df[["chrom", "start", "end", "clone", "cn_a", "cn_b"]].drop_duplicates()

    df["clone"] = df["clone"].astype(int)

    df = df[df["clone"].isin(clones)]

    df.to_csv(args.out_file, compression="gzip", index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--clone-file", required=True)

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
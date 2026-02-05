import pandas as pd


def main(args):
    
    df = pd.read_csv(args.tumour_content_file, sep="\t")

    df.add_prefix("tumour_content_", axis=1)

    df_e = pd.read_csv(args.evidence_file, index_col="run_type", sep="\t")

    df["normal_evidence"] = df_e.loc["normal", "evidence"]

    df["full_evidence"] = df_e.loc["full", "evidence"]

    df["bayes_factor"] = df["full_evidence"] - df["normal_evidence"]

    df.insert(0, 'patient', args.patient)
    
    df.insert(1, 'sample', args.sample)
    
    df.insert(2, 'coverage', args.coverage)

    df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--evidence-file", required=True)
    
    parser.add_argument("-t", "--tumour-content-file", required=True)
    
    parser.add_argument("-o", "--out-file", required=True)
    
    parser.add_argument("--patient", required=True, type=str)
    
    parser.add_argument("--sample", required=True, type=str)
    
    parser.add_argument("--coverage", required=True, type=float)
    
    cli_args = parser.parse_args()
    
    main(cli_args)

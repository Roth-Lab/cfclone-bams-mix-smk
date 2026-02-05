import pandas as pd


def df_to_wig(df: pd.DataFrame, genome_build: str, output_path: str) -> None:
    """
    Converts dataframe to a WIG fixedStep format file. Assumes input df is 
    hg38 build.

    Args:
        df (pd.DataFrame): DataFrame with columns ['chrom', 'start', 'end', 'reads']
        output_path (str): Path to write the .wig file
        
    Returns:
        None
    """
    with open(output_path, "w") as f:
        for c in df['chrom'].unique():
            df_c = (
                df
                .loc[df['chrom'] == c]
                .sort_values("start")
            )
            
            start = df_c["start"].iloc[0]
            end = df_c['end'].iloc[0]
            step = int(end - start)
            
            if start == 0:
                start = 1
            
            if genome_build == 'hg19':
                chrom = c.replace('chr')
            else:
                chrom = c
            
            s = "fixedStep chrom={chrom} start={start} step={step} span={step}\n".format(
                start=start,
                step=step,
                chrom=chrom
            )
            f.write(s)
            
            for read in df_c['reads'].values:
                f.write("{reads}\n".format(reads=int(read)))


def main(args):
    df = pd.read_csv(args.in_file, sep="\t")
    
    df_to_wig(df, genome_build='hg38', output_path=args.out_file)



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)

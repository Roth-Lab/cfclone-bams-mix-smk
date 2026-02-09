import pandas as pd

import pysam


def main(args):
    
    num_reads = pysam.view("-c", args.in_file).strip()
    
    try:
         
        num_reads = int(num_reads)
        
        obs_cov = (args.read_length * num_reads) / args.genome_length
    
    except:
        
        num_reads = 'Nan'
        
        obs_cov = "Nan"
    
    df = pd.DataFrame(data={'actual_coverage': [obs_cov]})
    
    df.insert(0, 'patient', args.patient)
    
    df.insert(1, 'sample', args.sample)
    
    df.insert(2, 'bam_id', args.bam_id)
    
    df.insert(3, 'coverage_id', args.coverage_id)
    
    df.insert(4, 'proportion_id', args.proportion_id)
    
    df.insert(5, 'coverage', args.coverage)
    
    df.insert(6, 'proportion', args.proportion)
    
    df.to_csv(args.out_file, sep='\t', index=False)
    

if __name__ == "__main__":
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    
    parser.add_argument('-i', '--in-file', type=str, required=True)
    
    parser.add_argument('-o', '--out-file', type=str, required=True)
    
    parser.add_argument('--patient', type=str, required=True)
    
    parser.add_argument('--sample', type=str, required=True)
    
    parser.add_argument('--coverage', type=float, required=True)
    
    parser.add_argument("--coverage-id", type=int, required=True)
    
    parser.add_argument("--proportion", type=float, required=True)
    
    parser.add_argument("--proportion-id", type=int, required=True)
    
    parser.add_argument("--bam-id", type=str, required=True)

    parser.add_argument('--read-length', type=int, required=True)
    
    parser.add_argument('--genome-length', type=int, required=True)
    
    cli_args = parser.parse_args()
    
    main(cli_args)
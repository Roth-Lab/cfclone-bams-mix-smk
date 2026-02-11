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
    
    samples = df_plot['final_sample'].unique()
    
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
    
    if df_patient_tc is not None:
        
        for s in ['initial_sample', 'final_sample']:
            
            patient = df_tmp['patient'].values[0]
            
            sample = df_tmp[s].values[0]
            
            df = df_patient_tc.loc[
                (df_patient_tc['patient_id'] == patient) & (df_patient_tc['sample_id'] == sample),
                ['mean', 'lower_hdi', 'upper_hdi']
            ]
            
            if s == 'initial_sample':
                
                x = 0.
            
            else:
                
                x = 1.
                
            ax.errorbar(
                x=x,
                y=df['mean'],
                yerr=[
                    df['mean'] - df['lower_hdi'],
                    df['upper_hdi'] - df['mean'],
                ],
                fmt="o",
                c='orange',
            )
    
    
    sample = df_tmp['final_sample'].values[0]
    
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
    
    patients = df['patient'].unique()
    
    initial_samples = df['initial_sample'].unique()
    
    final_samples = df['final_sample'].unique()
    
    assert len(patients) == 1
    
    assert len(initial_samples) == 1
    
    assert len(final_samples) == 1
    
    if args.cohort_tc_summary_file is not None:
        
        df_cohort = pd.read_csv(args.cohort_tc_summary_file, sep='\t')
        
        patients_cohort = df_cohort['patient_id'].unique()
        
        samples_cohort = df_cohort['sample_id'].unique()
        
        assert set(list(patients)) < set(list(patients_cohort))
        
        samples = list(initial_samples) + list(final_samples)
        
        assert set(samples) < set(list(samples_cohort))
        
        df_cohort = (
            
            df_cohort
            
            .loc[lambda df: (df['patient_id'] == patients[0]) & (df['sample_id'].isin(samples))]
            
            .drop_duplicates(subset=['patient_id', 'sample_id'])  # just use the first run
        )
        
    else:
        
        df_cohort = None
        
    plot(
        df_plot=df,
        df_patient_tc=df_cohort,
        output_path=args.out_file,
    )
    

if __name__ == "__main__":
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    
    parser.add_argument('-i', '--in-file', type=str, required=True) 
    
    parser.add_argument('-o', '--out-file', type=str, required=True)
    
    parser.add_argument("--cohort-tc-summary-file", type=lambda x: None if x == 'None' else x, required=True)

    cli_args = parser.parse_args()
    
    main(cli_args)
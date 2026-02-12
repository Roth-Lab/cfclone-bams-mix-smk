import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt


plot_settings = {
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "figure.titleweight": "bold",
    "figure.titlesize": 12,
}


def plot(
    df_plot: pd.DataFrame,
    df_patient_tc: pd.DataFrame | None = None,
    figsize: tuple[int, int] | None = None,
    output_path: str | None = None,
    show: bool = False
) -> None:
    with plt.rc_context(plot_settings):
    
        patients = df_plot['patient'].unique()
        
        num_patients = len(patients)
        
        assert num_patients == 1
        
        initial_samples = df_plot['initial_sample'].unique()
        
        num_initial_samples = len(initial_samples)
        
        assert num_initial_samples == 1
    
        coverage_ids = df_plot['coverage'].unique()
        
        num_coverage_ids = len(coverage_ids)
        
        assert num_coverage_ids == 1
        
        samples = df_plot['final_sample'].unique()
        
        num_samples = len(samples)
    
        num_cols = int(np.ceil(np.sqrt(num_samples)))
        
        num_rows = int(np.ceil(np.sqrt(num_samples / num_cols)))
        
        figsize = (num_cols * 6, num_rows * 4)
            
        fig = plt.figure(figsize=figsize, constrained_layout=True)
        
        gs = fig.add_gridspec(num_rows, num_cols)
        
        for row_idx in range(num_rows):
            
            for col_idx in range(num_cols):
                
                ax = fig.add_subplot(gs[row_idx, col_idx])
                
                sample_idx = row_idx * num_cols + col_idx
                
                if sample_idx < num_samples:
                
                    df_sample = df_plot.loc[df_plot['final_sample'] == samples[sample_idx]]
                    
                    plot_sample(
                        df_sample=df_sample,
                        ax=ax,
                        df_patient_tc=df_patient_tc,
                    )
                    
                else:
                    
                    ax.axis('off')
                    
        fig.suptitle("Patient: {p}\n Initial sample:{i}\n Coverage: {c}".format(p=patients[0], i=initial_samples[0], c=coverage_ids[0]))
                
        if output_path:
            
            plt.savefig(output_path)

        if show:
            
            plt.show()

        plt.close()
        
    
def plot_sample(
    df_sample: pd.DataFrame,
    ax: plt.Axes,
    df_patient_tc: pd.DataFrame | None = None,
):
    df_sample = df_sample.sort_values(by=['proportion'])
    
    ax.scatter(
        x=df_sample['proportion'], 
        y=df_sample["mean"]
    )
    
    ax.errorbar(
        x=df_sample['proportion'],
        y=df_sample["mean"],
        yerr=[
            df_sample["mean"] - df_sample["lower_hdi"],
            df_sample["upper_hdi"] - df_sample["mean"],
        ],
        fmt="o",
        ecolor="gray",
        alpha=0.5,
    )
    
    if df_patient_tc is not None:
        
        for s in ['initial_sample', 'final_sample']:
            
            patient = df_sample['patient'].values[0]
            
            sample = df_sample[s].values[0]
            
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
    
    ax.set_ylabel("Tumour Content estimate")
    
    ax.set_xlabel("Proportion of final sample".format(sample))
    
    ax.set_title("Final sample: {}".format(sample))

    
def main(args):
    
    df = pd.read_csv(args.in_file, sep='\t')
    
    patients = df['patient'].unique()
    
    initial_samples = df['initial_sample'].unique()
    
    final_samples = df['final_sample'].unique()
    
    assert len(patients) == 1
    
    assert len(initial_samples) == 1
    
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
    
    parser.add_argument('-i', '--in-file', type=str, default='summary_new_new.tsv') 
    
    parser.add_argument('-o', '--out-file', type=str, default='lol.png')
    
    parser.add_argument("--cohort-tc-summary-file", type=lambda x: None if x == 'None' else x, default='/home/matteo/projects/cfdna/data/cohort/tumour_contents.tsv')

    cli_args = parser.parse_args()
    
    main(cli_args)
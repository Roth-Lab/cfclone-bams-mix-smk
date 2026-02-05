import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

colour_blind_friendly = {
    'blue': '#377eb8',
    'orange': '#ff7f00',
    'green': '#4daf4a',
    'pink': '#f781bf',
    'brown': '#a65628',
    'purple': '#984ea3',
    'grey': '#999999',
    'red': '#e41a1c',
    'yellow': '#dede00'
}

plot_settings = {
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.labelweight": "bold",
}


def plot(
    df_plot: pd.DataFrame,
    df_patient_bf: pd.DataFrame | None = None,
    xlims: tuple[float, float] | None = None,
    figsize: tuple[int, int] | None = None,
    output_path: str | None = None,
    show: bool = False
) -> None:
    with plt.rc_context(plot_settings):
        _plot(
            df_plot=df_plot,
            df_patient_bf=df_patient_bf,
            xlims=xlims,
            figsize=figsize,
            output_path=output_path,
            show=show
        )
    
def _plot(
    df_plot: pd.DataFrame,
    df_patient_bf: pd.DataFrame | None = None,
    xlims: tuple[float, float] | None = None,
    figsize: tuple[int, int] | None = None,
    output_path: str | None = None,
    show: bool = False
) -> None:
    
    sample_ids = df_plot['sample'].unique()
    
    num_sample_ids = len(sample_ids)
    
    num_cols = int(np.ceil(np.sqrt(num_sample_ids)))
    
    num_rows = int(np.ceil(num_sample_ids / num_cols))
    
    if figsize is None:
        
        figsize = (num_cols * 6, num_rows * 4)
        
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    
    gs = fig.add_gridspec(num_rows, num_cols)
    
    for i in range(num_rows):
        
        for j in range(num_cols):
            
            index = i * num_cols + j
            
            if index >= num_sample_ids:
                
                continue
            
            sample_id = sample_ids[index]
            
            df_tmp = df_plot.loc[df_plot['sample'] == sample_id]
            
            ax = fig.add_subplot(gs[i, j])
            
            c = get_colours(df_tmp['bayes_factor'])
            
            # if sample_id == 'VOA10055P':
                
            ax.scatter(x=df_tmp["coverage"], y=df_tmp["bayes_factor"], c=c)
            
            if df_patient_bf is not None:
                
                df = df_patient_bf.loc[df_patient_bf['sample'] == sample_id]
                
                c = get_colours(df['bayes_factor'])
                
                ax.scatter(x=df['coverage'], y=df['bayes_factor'], c=c, marker='x')
            
            ax.set_ylabel("Bayes Factor (log)")
            
            ax.set_xlabel("Coverage")
            
            if xlims is not None:
                
                ax.set_xlim(xlims)
            
            ax.xaxis.set_inverted(True)
            
            ax.set_xscale('log')
            
            ax.set_title("Sample {}".format(sample_id))
            
            handles, labels = add_legend()
            
            ax.legend(handles=handles, labels=labels, loc='best')
    
    if output_path:
        
        plt.savefig(output_path)

    if show:
        
        plt.show()   
        
    plt.close()

 
def add_legend() -> tuple[list, list]:
    handles = []
    labels = []
    
    handles.append(
        plt.Line2D(
            xdata=[0],
            ydata=[0],
            marker='o',
            color='w',
            label='Detected',
            markerfacecolor=colour_blind_friendly['green'], 
            # markersize=8
        )
    )
    labels.append('Detected')
    
    handles.append(
        plt.Line2D(
            [0], [0], marker='o', color='w', label='Inconclusive',
            markerfacecolor=colour_blind_friendly['blue']
        )
    )
    labels.append('Inconclusive')
    
    handles.append(
        plt.Line2D(
            [0], [0], marker='o', color='w', label='Not Detected',
            markerfacecolor=colour_blind_friendly['red'],
        )
    )
    labels.append('Not Detected')
    
    return handles, labels
    
    
def get_colours(bayes_factors: np.ndarray) -> list[str]:
    c = []
    for bf in bayes_factors:
        if bf > 3.:
            # detected tumour fragment 
            c.append(colour_blind_friendly['green'])
        elif -3 < bf and bf < 3:
            # Inconclusive
            c.append(colour_blind_friendly['blue'])
        else:
            # Not detected
            c.append(colour_blind_friendly['red'])
    return c


def main(args):
    
    df = pd.read_csv(args.in_file, sep='\t')
    
    plot(df_plot=df, output_path=args.out_file)



if __name__ == "__main__":
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    
    parser.add_argument("-i", "--in-file", required=True)
    
    parser.add_argument("-o", "--out-file", required=True)
    
    cli_args = parser.parse_args()
    
    main(cli_args)
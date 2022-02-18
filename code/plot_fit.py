
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
sn.set()

from compute_5PL import r_squared, PL5, fit_curve, get_standard


def plot_data(ax, s, conc, fit, color="red"):
    _ = ax.scatter(conc, fit, color=color, s=.5, alpha=0.5)
    _ = ax.scatter(s['conc'], s['FI'], color=color, s=45)



def plot_s(data_df, col_df, gene="M-CSF", run="20211021", fit=True, params=[]):
    '''
    retrieve the standard curve and plot along with fit line
    '''
    # init the plot
    fig, ax = plt.subplots(figsize=(12,12))
    
    
    s = get_standard(data_df, col_df, run=run, gene=gene)
    
    
    # fit the curve if no params are given
    if fit:
        if len(params):
            # use params to plot the curve
            R = r_squared(params, s)
            conc = np.power(10, np.linspace(-3, 5.1, 10000))
            fit = PL5(conc, params)  
        else:
            # fit the curve
            conc, fit, R, params = fit_curve(s)
            params = list(params)
        
    plot_data(ax, s, conc, fit, color="blue")
    
    # plots the params info
    ymax = s['FI'].max()
    x = 0.001
    for i,p in enumerate("ABCDE"):
        _ = ax.text(x,ymax - i*ymax/10, f"{p}: {round(params[i],3)}", size=30, color="darkgray")
    
    _ = ax.set_ylim(0,ymax*1.1)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=20)
    _ = ax.set_ylabel('FI', fontsize=20)
    _ = plt.xticks(fontsize=20)
    _ = plt.yticks(fontsize=20)
    plt.title(f"{gene} | R={round(R,5)}", fontsize=30)
    
    return


def plot_multi_s(data_df, col_df, gene="M-CSF", runs=["20211021", "20211222"]):
    
    fig, ax = plt.subplots(figsize=(12,12))
    
    c = {"20211021":"red", "20211222":"green"}
    
    ymax = 0
    for r in runs:
        s = get_standard(data_df, col_df, run=r, gene=gene)
        # fit the curve
        conc, fit, R, _ = fit_curve(s)
        plot_data(ax, s, conc, fit, color=c[r])
        ymax = max(ymax, s['FI'].max())
    _ = ax.set_ylim(0,ymax*1.1)
    _ = ax.set_xscale('log')
    _ = ax.set_xlabel('conc', fontsize=20)
    _ = ax.set_ylabel('FI', fontsize=20)
    _ = plt.xticks(fontsize=20)
    _ = plt.yticks(fontsize=20)
    _ = ax.set_xticks([10**i for i in range(-3,6)])
    plt.title(f"{gene} | {'vs '.join(runs)}", fontsize=30)
    
    return

import random
import entropy_estimators as ee
from scipy import stats
from os import listdir
import matplotlib.pyplot as plt
import scipy.io as sio
import matplotlib.cm as cm
import numpy as np
import matplotlib.patches as mpatches
import scipy.stats as sp
import pandas as pd
from scipy.optimize import curve_fit
import seaborn as sns
import matplotlib

def variance_t(t,tau_x,a):
    return a*(tau_x/2)*(1-np.exp(-t/tau_x)**2)
def var_rate(t, g, tx, ty,z):
    a = (1 / ((tx - ty) ** 2 * (tx + ty))) * np.exp(-t * (1 / tx + 1 / ty))
    b = 4 * g ** 2 * tx ** 2 * ty ** 3
    c = np.exp(t * (-(1 / tx) + 1 / ty)) * g ** 2 * tx ** 2 * ty ** 2 * (tx + ty)
    d = np.exp(t * (1 / tx + 1 / ty)) * (tx - ty) ** 2 * (tx + ty + g ** 2 * tx * ty ** 2)
    e = np.exp(t * (1 / tx - 1 / ty)) * (tx + ty) * (tx ** 2 + ty ** 2 + tx * ty * (-2 + g ** 2 * ty ** 2))

    return z*a * (b - c + d - e)

if __name__ == "__main__":
    #white backgrounds
    sns.set_style('white')
    matplotlib.rcdefaults()
    matplotlib.rcParams['pdf.fonttype'] = 42
    # matplotlib.rcParams['svg.fonttype'] = 42

    ### load CRISPR data -- Wild type then CRISPRi strain (14th bp mutation)
    folder1 = 'Fig2_SourceData_1'
    names1 = listdir(folder1)

    if '.DS_Store' in names1:
        names1.remove('.DS_Store')
    names1 = sorted(names1, key=lambda x: int(x.split('.')[0]))
    names1 = [folder1 + '/' + name for name in names1]


    folder2 = 'Fig2_SourceData_2'
    names2 = listdir(folder2)

    if '.DS_Store' in names2:
        names2.remove('.DS_Store')
    names2 = sorted(names2, key=lambda x: int(x.split('.')[0]))
    names2 = [folder2 + '/' + name for name in names2]

    names=names1+names2

    #initialize dataframe for presenting data
    df = pd.DataFrame({'X': [], 'Y': [], 'Z': [], 'Information': [], 'Correlation': [], 'nCells': [], 'label': []})
    tags = ['WT', 'CRISPRi Lon']
    n = 0

    #set number of samples
    num_samps_each = 5

    total_x_wt = np.array([])
    total_x_c = np.array([])

    total_y_wt = np.array([])
    total_z_wt = np.array([])
    total_y_c = np.array([])
    total_z_c = np.array([])

    for j, name in enumerate(names):
        print(name)
        mat_contents = sio.loadmat(name)
        color_mat = mat_contents['data3D']
        t_len = mat_contents['data3D'].shape[2]

        for t in range(t_len):
            x = color_mat[:, 5, t][np.where(~np.isnan(color_mat[:, 5, t]))]
            y = color_mat[:, 7, t][np.where(~np.isnan(color_mat[:, 5, t]))]
            z = color_mat[:, 9, t][np.where(~np.isnan(color_mat[:, 5, t]))]
            spearman = sp.spearmanr(color_mat[:, 7, t][np.where(~np.isnan(color_mat[:, 5, t]))],
                                    color_mat[:, 9, t][np.where(~np.isnan(color_mat[:, 5, t]))])[0]
            if n < 2:
                if n == 0:
                    if t == t_len - 1:
                        total_x_wt = np.append(total_x_wt, x.ravel())
                        total_y_wt = np.append(total_y_wt, y.ravel())
                        total_z_wt = np.append(total_z_wt, z.ravel())
                if n == 1:
                    if t == t_len - 1:
                        total_x_c = np.append(total_x_c, x.ravel())
                        total_y_c = np.append(total_y_c, y.ravel())
                        total_z_c = np.append(total_z_c, z.ravel())
            try:
                info = ee.mi(ee.vectorize(y), ee.vectorize(z), k=3)
            except:
                info = 0

            n_cells = len(np.where(~np.isnan(color_mat[:, 5, t]))[0])
            df2 = pd.DataFrame(
                {'X': [np.std(x) / np.mean(x)], 'Y': [np.std(y) / np.mean(y)], 'Z': [np.std(z) / np.mean(z)],
                 'Information': [info], 'Correlation': [spearman], 'nCells': [n_cells], 'label': [tags[n]]})
            frames = [df, df2]
            df = pd.concat(frames)
        if (j + 1) % num_samps_each == 0:
            n += 1


    ### Scatter plots
    fnaut, axnaut = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(2, 5))
    axnaut[0].set_xlim([400, 2500])
    axnaut[0].set_ylim([400, 1200])
    axnaut[0].set_xlabel('pinaA')
    axnaut[1].set_xlabel('pinaA')
    axnaut[1].set_ylabel('pacrAB')
    axnaut[0].set_ylabel('pacrAB')

    axnaut[0].set_title('Crispr Lon')
    axnaut[1].set_title('Wild Type')
    colors = ['darkblue', 'salmon']
    axnaut[1].scatter(total_y_wt, total_z_wt, color=colors[0], alpha=0.1)
    axnaut[0].scatter(total_y_c, total_z_c, color=colors[1], alpha=0.1)
    fnaut.subplots_adjust(hspace=.5)
    fnaut.savefig('fig_out/crisp_2dhist.pdf', bbox_inches='tight')

    ### Histograms
    plt.close('all')
    fnaut2, axnaut2 = plt.subplots(figsize=(2, 5))
    sns.distplot(total_x_wt, ax=axnaut2, color='darkblue', hist=False, kde_kws={"lw": 3}, vertical=True)
    sns.distplot(total_x_c, ax=axnaut2, color='salmon', hist=False, kde_kws={"lw": 3}, vertical=True)
    stats.ks_2samp(total_x_wt, total_x_c)
    axnaut2.set_xlabel('Probability density')
    axnaut2.set_ylabel('Fluorescence (A.U.)')
    fnaut2
    fnaut2.savefig('fig_out/crisp_hist.pdf', bbox_inches='tight')


    ### Time series
    plt.close('all')

    f, ax = plt.subplots(1, 2, sharex=True, figsize=(10, 3))
    fnaut, axnaut = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(4, 10))
    axnaut[0].set_xlim([400, 5000])
    axnaut[0].set_ylim([400, 5000])
    ax[0].set_xlim([-10, 300])
    sns.set_style('white')
    colors = ['darkblue', 'salmon']
    tags = ['WT', 'CRISPRi Lon']
    tvect = np.linspace(0, 400, 1000)
    recs = []

    for j, color in enumerate(colors):
        x, y = df.loc[df.label == tags[j]].nCells, df.loc[df.label == tags[j]].X
        x, y2 = df.loc[df.label == tags[j]].nCells, df.loc[df.label == tags[j]].Information
        ys = df.loc[df.label == tags[j]].Y
        zs = df.loc[df.label == tags[j]].Z

        nbins = 15
        bins = np.linspace(0, 275, nbins)
        idx = np.digitize(x, bins)
        means = []
        errors = []
        means2 = []
        errors2 = []
        for i in range(nbins):
            ax[0].errorbar(bins[i] + (200 / (2 * nbins)), np.mean(y[idx == i + 1]), fmt='o', color=color)

            ax[1].errorbar(bins[i] + (200 / (2 * nbins)), np.mean(y2[idx == i + 1]), fmt='o', color=color)

            #           ax.errorbar(bins[i]+(200/(2*nbins)),np.mean(y[idx==i+1]),yerr=np.std(y[idx==i+1]),fmt='o',color=color)
            means.append(np.mean(y[idx == i + 1]))
            errors.append(sp.sem(y[idx == i + 1]))
            means2.append(np.mean(y2[idx == i + 1]))
            errors2.append(sp.sem(y2[idx == i + 1]))

        means = np.asarray(means)
        errors = np.asarray(errors)
        ax[0].fill_between(bins + (200 / (2 * nbins)), means - errors, means + errors, color=color, alpha=0.5)
        means2 = np.asarray(means2)
        errors2 = np.asarray(errors2)
        ax[1].fill_between(bins + (200 / (2 * nbins)), means2 - errors2, means2 + errors2, color=color, alpha=0.5)
        recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=colors[j]))

    ax[0].legend(recs, ['Wild Type', 'CRISPR Lon'], title='Strain', loc=4)
    ax[1].legend(recs, ['Wild Type', 'CRISPR Lon'], title='Strain', loc=4)

    ax[0].set_ylabel('Coefficient of Variation')
    ax[0].set_title('Activator Variance over Time')
    ax[0].set_xlabel('Number of Cells in Microcolony')
    ax[1].set_xlabel('Number of Cells in Microcolony')

    ax[1].set_title('Downstream coordination over time')
    ax[1].set_ylabel('Infromation (bits)')

    plt.tight_layout()

    f.savefig('fig_out/crips_ts.png', bbox_inches='tight')

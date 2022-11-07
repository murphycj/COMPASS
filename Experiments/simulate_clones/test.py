# %%
import glob
from pathlib import Path

import pandas as pd
from scipy import stats
import numpy as np
from scipy.cluster import hierarchy
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %% [markdown]
# # variant test

# %%
def variant_prob(obs, N, mixing, background):
    prob_het = stats.binom.pmf(obs, N, 0.5)
    prob_hom = stats.binom.pmf(obs, N, 0.95)
    prob_variant = 1 - (prob_het * prob_hom) - (1 - prob_het) * (1 - prob_hom)
    
    prob_mix_1 = stats.binom.pmf(obs, N, 0.33)
    prob_mix_2 = stats.binom.pmf(obs, N, 0.66)
    prob_mix = 1 - (prob_mix_1 * prob_mix_2) - (1 - prob_mix_1) * (1 - prob_mix_2)
    
    real = (1 - mixing) * prob_variant + mixing * prob_mix
    noise = stats.binom.pmf(obs, N, background)

    llr = np.log(real) - np.log(noise)
    
    return pd.DataFrame([{
        'obs': obs,
        'N': N,
        'mixing': mixing,
        'background': background,
        'prob_variant': prob_variant,
        'prob_mix': prob_mix,
        'noise': noise,
        'prob_real': real,
        'prob_noise': noise,
        'llr': llr
    }])

# %%
probs = pd.concat([variant_prob(i+1, 100, 0.1, 0.1) for i in range(100)])

# %%
for col in ['llr', 'prob_variant', 'prob_mix', 'noise', 'prob_real', 'prob_noise', 'llr']:
    fig = px.scatter(data_frame=probs, x='obs', y=col)
    fig.show(renderer='png')


# %%
N = 50
obs = 4
mixing = 0.1
background = 0.01

real = (1 - mixing) * stats.binom.pmf(obs, N, 0.5) + mixing * stats.binom.pmf(obs, N, 0.33)
noise = stats.binom.pmf(obs, N, background)

llr = np.log(real) - np.log(noise)

real, noise, llra
astats.binom.rvs(10, 0.99)

# %%
np.log(real)

# %%
background = 0.001
cells = 3000
obs = 30

pval = 1 - stats.binom.cdf(obs, cells, background)a

pval, np.log(pval)

# %%
stats.betabinom.rvs(10, 0.1, 1, size=10)

# %%
data = stats.binom.rvs(cells, [0.001, 0.01], size=(5, 2))
data

# %% [markdown]
# # compass

# %%
def parse_af(x):
    ref, alt, _ = x.split(':')
    ref = int(ref)
    alt = int(alt)
    total = ref + alt
    
    if total == 0:
        return 0
    else:
        return alt / total
    

def parse_input(fpath):
    test = pd.read_csv(fpath)
    
    idx = np.where(test.columns == '0')[0][0]
    
    test = test[test.columns[idx:]]
    print(test)
    test = test.applymap(lambda x: parse_af(x))
    test.columns = test.columns.map(lambda x: f'Cell {x}')
    
    return test
    

def parse_output(cell_assignment, genotypes):
    cells = pd.read_csv(cell_assignment, sep='\t')
    cells['cell'] = cells['cell'].map(lambda x: f'Cell {x}')
    genotypes = pd.read_csv(genotypes, sep='\t', index_col=0)

    mut = pd.DataFrame(3, index=genotypes.columns, columns=cells['cell'])

    for i in cells.index:
        cell = cells.at[i, 'cell']
        node = cells.at[i, 'node']
        node = f'Node {node}'
        
        for j in genotypes.columns:
            try:
                mut.at[j, cell] = int(genotypes.at[node, j])
            except:
                pass
    
    return mut

def comparison(mut1_fpath, mut2_ca_fpath, mut2_g_fpath):
    mut1 = parse_input(mut1_fpath)

    mut2 = parse_output(mut2_ca_fpath, mut2_g_fpath)

    row_link = hierarchy.linkage(mut2.T, method="ward")
    row_dendro = hierarchy.dendrogram(row_link, no_plot=True)
    row_order = row_dendro["leaves"][::-1]

    mut2 = mut2[mut2.columns[row_order]]
    mut1 = mut1[mut1.columns[row_order]]

    freq1 = round(100 * (((mut1>=0.2).sum(axis=0)==2).sum() / mut1.shape[1]), 2)
    freq2 = round(100 * (((mut2.isin([1,2])).sum(axis=0)==2).sum() / mut2.shape[1]), 2)
    
    return mut1, mut2, freq1, freq2


def plot(fout, mut1, mut2, freq1, freq2):
    fig = make_subplots(rows=2, cols=1, row_titles=['Raw (VAF)', 'COMPASS (genotypes)'])

    fig.add_trace(
        go.Heatmap(
            z=mut1*100,
            y=mut2.index,
            colorbar=dict(y=.8,len=.5)
        ),
        row=1,
        col=1
    )

    fig.add_trace(
        go.Heatmap(
            z=mut2,
            y=mut2.index,
            colorbar=dict(y=.2,len=.5)
        ),
        row=2,
        col=1
    )

    # fig = px.imshow(mut2, aspect='auto')
    fig.update_layout(title=f'Expected clone frequency: {freq1}%<br>Observed clone frequency: {freq2}%')
    fig.update_xaxes(showticklabels=False, title='Cells')
    fig.update_yaxes(title='Mutations')
    # fig.show(renderer='png', width=1200, height=600)
    fig.write_image(fout, width=1200, height=600)
    
    
def plot_single(mut_ca_fpath, mut_g_fpath):
    mut = parse_output(mut_ca_fpath, mut_g_fpath)
    
    row_link = hierarchy.linkage(mut.T, method="ward")
    row_dendro = hierarchy.dendrogram(row_link, no_plot=True)
    row_order = row_dendro["leaves"][::-1]

    mut = mut[mut.columns[row_order]]
    
    fig = px.imshow(mut, aspect='auto')
    fig.update_layout(title=f'Expected clone frequency: {freq1}%<br>Observed clone frequency: {freq2}%')
    fig.update_xaxes(showticklabels=False, title='Cells')
    fig.update_yaxes(title='Mutations')
    fig.show(renderer='png', width=1200, height=600)
    # fig.write_image(fout, width=1200, height=600)

# %%
plot_single(
    '/Users/charliemurphy/Desktop/github/forks/COMPASS/data/AML-99-001_cellAssignments.tsv',
    '/Users/charliemurphy/Desktop/github/forks/COMPASS/data/AML-99-001_nodes_genotypes.tsv'
)

# %%
np.where(pd.read_csv('/Users/charliemurphy/Desktop/github/forks/COMPASS/data/AML-99-001_variants.csv').columns == '0')[0][0]

# %%

# mut1, mut2, freq1, freq2 = comparison(
#     '/Users/charliemurphy/Desktop/github/forks/COMPASS/Experiments/simulate_clones/test_variants.csv',
#     '/Users/charliemurphy/Desktop/github/forks/COMPASS/Experiments/simulate_clones/test_cellAssignments.tsv',
#     '/Users/charliemurphy/Desktop/github/forks/COMPASS/Experiments/simulate_clones/test_nodes_genotypes.tsv'
# )


mut1, mut2, freq1, freq2 = comparison(
    '/Users/charliemurphy/Desktop/github/forks/COMPASS/data/AML-99-001_variants.csv',
    '/Users/charliemurphy/Desktop/github/forks/COMPASS/data/AML-99-001_cellAssignments.tsv',
    '/Users/charliemurphy/Desktop/github/forks/COMPASS/data/AML-99-001_nodes_genotypes.tsv'
)

fig = make_subplots(rows=2, cols=1, row_titles=['Raw (VAF)', 'COMPASS (genotypes)'])

fig.add_trace(
    go.Heatmap(
        z=mut1*100,
        y=mut2.index,
        colorbar=dict(y=.8,len=.5)
    ),
    row=1,
    col=1
)

fig.add_trace(
    go.Heatmap(
        z=mut2,
        y=mut2.index,
        colorbar=dict(y=.2,len=.5)
    ),
    row=2,
    col=1
)

# fig = px.imshow(mut2, aspect='auto')
fig.update_layout(title=f'Expected clone frequency: {freq1}%<br>Observed clone frequency: {freq2}%')
fig.update_xaxes(showticklabels=False, title='Cells')
fig.update_yaxes(title='Mutations')
fig.show(renderer='png', width=1200, height=600)

# %%
basedir = Path('/Users/charliemurphy/Desktop/github/forks/COMPASS/Experiments/simulate_clones/')

res = []

for fpath in glob.glob(str(basedir / 'COMPASS_out/*_nodes_genotypes.tsv')):
    name = Path(fpath).stem.replace('_nodes_genotypes', '')
    nnodes, freq, rep = name.split('_')
    nnodes = int(nnodes)
    freq = [float(i) for i in freq.split('-')]
    
    mut1, mut2, freq1, freq2 = comparison(
        basedir / f'generate/{name}/{name}_variants.csv',
        basedir / f'COMPASS_out/{name}_cellAssignments.tsv',
        basedir / f'COMPASS_out/{name}_nodes_genotypes.tsv'
    )
    
    plot(basedir / f'COMPASS_out/{name}_heatmap.png', mut1, mut2, freq1, freq2)
    
    res.append({
        'name': name,
        'nnodes': nnodes,
        'freq': ';'.join([str(i * 100) for i in freq]),
        'wt_freq': freq[0] * 100,
        'mut_freq': 100 - (freq[0] * 100),
        'freq_exp': freq1,
        'freq_obs': freq2,
        'freq_diff': freq2 - freq1
    })
    
res = pd.DataFrame(res)
res['nnodes'] = res['nnodes'].astype(str)

# %%


# %%
fig = px.scatter(data_frame=res, x='freq_exp', y='freq_obs', color='nnodes')
fig.update_yaxes(type='log', title='Clone frequency, expected (%)')
fig.update_xaxes(type='log', title='Clone frequency, observed (%)')
fig.show(renderer='png', width=1000)

fig = px.box(data_frame=res, x='mut_freq', y='freq_diff', color='nnodes', points='all')
fig.show(renderer='png', width=1000)

fig = px.box(data_frame=res, x='mut_freq', y='freq_obs', color='nnodes', points='all')
fig.show(renderer='png', width=1000)


# %%
res

# %%




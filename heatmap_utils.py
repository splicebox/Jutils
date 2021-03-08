import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_heatmap(file, meta_file, out_dir, p_value_threhold, q_value_threhold,
                 dpsi_threshold, foldchange_threshold, avg_threshold,
                 unsupervised, aggregate, method, metric, prefix, top):

    samples = []
    conditions = []
    sample_cond_dict = {}
    with open(meta_file, 'r') as f:
        for line in f:
            sample, cond = line.strip().split('\t')
            samples.append(sample)
            conditions.append(cond)
            sample_cond_dict[sample] = cond

    data_df = pd.read_csv(file, sep='\t', comment='#')
    data_df['GeneName'] = data_df['GeneName'].apply(truncate_gene_name)
    original_columns = data_df.columns
    if np.all(data_df['ReadCount1'] == '.') and np.all(data_df['PSI'] == '.'):
        raise Exception("Column 'ReadCount1' and 'PSI' don't contain values! This could happen if the TSV file is generated from MAJIQ outputs.")
        sys.exit('-1')

    if unsupervised:
        data_df = process_data_unsupervised(data_df, samples, original_columns, avg_threshold, aggregate, top)
    else:
        data_df = process_data_supervised(data_df, samples, sample_cond_dict, conditions, original_columns,
                                          p_value_threhold, q_value_threhold, dpsi_threshold, foldchange_threshold,
                                          avg_threshold, aggregate)

    generate_heatmaps(data_df, original_columns, conditions, sample_cond_dict,
                      method, metric, prefix, aggregate, unsupervised, out_dir)


def truncate_gene_name(string):
    gene_names = string.split(',')
    if len(gene_names) > 3:
        return ','.join(gene_names[:3])
    return string


def process_group(x):
    start, end = float('inf'), 0
    for label in x['FeatureLabel']:
        _chr, coord = label.split(':')
        _start, _end = (int(v) for v in coord.split('-'))
        start, end = min(start, _start), max(end, _end)
    return f"{x['GeneName'].iloc[0]}_{_chr}:{start}-{end}"


def process_data_unsupervised(data_df, samples, original_columns, avg_threshold, aggregate, top_num):
    print('Unsupervised mode, the option --q-value, --p-value, --dpsi, --fold-change will be ignored')
    if np.any(data_df['GeneName'] != '.'):
        data_df = data_df[data_df['GeneName'] != '.']

    # counts_df = data_df['ReadCount1'].str.split(',', expand=True).astype(float)
    # bools = counts_df.apply(lambda x: x.mean() > 5, axis=1)
    # data_df = data_df[bools]

    if 'dPSI' in original_columns:
        if data_df.shape[1] > 14 and aggregate:
            print("Multi-way comparison doesn't support aggregate mode, aggregate mode disabled!")
            aggregate = False

        if aggregate:
            data_df2 = data_df['PSI'].str.split(',', expand=True).astype(float)
            data_df2.columns = samples
            data_df2[['GeneName', 'FeatureLabel', 'GroupID', 'dPSI']] = data_df[['GeneName', 'FeatureLabel', 'GroupID', 'dPSI']]

            groups = data_df2.groupby(['GroupID'])
            data_df3 = groups.apply(lambda x: x.iloc[x[samples[0]].argmax()])
            data_df3['Bool'] = data_df3['dPSI'].apply(lambda x: 1 if x > 0 else -1)
            data_df4 = pd.merge(data_df2, data_df3['Bool'], on='GroupID')
            data_df4 = data_df4[data_df4['dPSI'] * data_df4['Bool'] > 0]

            groups = data_df4.groupby(['GroupID'])
            data_df5 = groups.apply(lambda x: np.sum(x[samples], axis=0))
            data_df5.index = groups.apply(process_group)
            # data_df5['variance'] = data_df5.mean(axis=1) / data_df5.var(axis=1)
            data_df5['variance'] = data_df5.sub(data_df5.mean(axis=1), axis=0).abs().mean(axis=1)
            data_df5 = data_df5.sort_values(by=['variance'], ascending=False)
            data_df = data_df5.drop(columns=['variance'])
        else:
            data_df2 = data_df['PSI'].str.split(',', expand=True).astype(float)
            data_df2.columns = samples

            # data_df2['variance'] = data_df2.mean(axis=1) / data_df2.var(axis=1)
            data_df2['variance'] = data_df2.sub(data_df2.mean(axis=1), axis=0).abs().mean(axis=1)
            data_df2['index'] = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)

            groups = data_df2.groupby(['index'])
            data_df = groups.apply(lambda x: x.iloc[x['variance'].argmax()])
            data_df = data_df.sort_values(by=['variance'], ascending=False)
            data_df = data_df.drop(columns=['variance', 'index'])

    elif 'log2FoldChange' in original_columns:
        if aggregate:
            print('Warming: the aggregate mode will be ignore for DSA data!')
        data_df = data_df[data_df.iloc[:, 12:].apply(lambda x: np.any(x > avg_threshold) ,axis=1)]
        data_df2 = data_df['ReadCount1'].str.split(',', expand=True).astype(float)
        data_df2.columns = samples
        data_df2['variance'] = data_df2.var(axis=1) / data_df2.mean(axis=1)
        data_df2 = data_df2.apply(lambda x: np.log2(x + 1e-2))
        # data_df2['variance'] = data_df2.mean(axis=1) / data_df2.var(axis=1)
        data_df2.index = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
        data_df = data_df2
        data_df = data_df.sort_values(by=['variance'], ascending=False)
        data_df = data_df.drop(columns=['variance'])

    num = min(top_num, data_df.shape[0])
    data_df = data_df.iloc[:num, :]
    return data_df


def process_data_supervised(data_df, samples, sample_cond_dict, conditions, original_columns,
                            p_value_threhold, q_value_threhold, dpsi_threshold,
                            foldchange_threshold, avg_threshold, aggregate):
    data_df = data_df[data_df['p-value'] < p_value_threhold]
    data_df = data_df[data_df['q-value'] < q_value_threhold]
    if np.any(data_df['GeneName'] != '.'):
        data_df = data_df[data_df['GeneName'] != '.']

    if 'dPSI' in original_columns:
        if data_df.shape[1] > 14 and aggregate:
            print("Multi-way comparison doesn't support aggregate mode, aggregate mode disabled!")
            aggregate = False

        if aggregate:
            selected_groups = data_df[abs(data_df['dPSI']) > dpsi_threshold]['GroupID'].drop_duplicates()
            data_df = data_df[data_df['GroupID'].isin(selected_groups)]

            data_df2 = data_df['PSI'].str.split(',', expand=True).astype(float)
            data_df2.columns = samples
            data_df2[['GeneName', 'FeatureLabel', 'GroupID', 'dPSI']] = data_df[['GeneName', 'FeatureLabel', 'GroupID', 'dPSI']]

            groups = data_df2.groupby(['GroupID'])
            cond1_samples = [sample for sample, cond in sample_cond_dict.items() if cond == conditions[0]]
            data_df3 = groups.apply(lambda x: x.iloc[x[cond1_samples].max(axis=1).argmax()])
            data_df3['Bool'] = data_df3['dPSI'].apply(lambda x: 1 if x > 0 else -1)
            data_df4 = pd.merge(data_df2, data_df3['Bool'], on='GroupID')
            data_df4 = data_df4[data_df4['dPSI'] * data_df4['Bool'] > 0]

            groups = data_df4.groupby(['GroupID'])
            data_df = groups.apply(lambda x: np.sum(x[samples], axis=0))
            data_df.index = groups.apply(process_group)
        else:
            data_df = data_df[abs(data_df['dPSI']) > dpsi_threshold]
            data_df2 = data_df['PSI'].str.split(',', expand=True).astype(float)
            data_df2.columns = samples
            data_df2['index'] = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
            data_df2['dPSI'] = np.abs(data_df['dPSI'])

            groups = data_df2.groupby(['index'])
            data_df = groups.apply(lambda x: x.iloc[x['dPSI'].argmax()])
            data_df = data_df.sort_values(by=['dPSI'], ascending=False)
            data_df = data_df.drop(columns=['dPSI', 'index'])

    elif 'log2FoldChange' in original_columns:
        if aggregate:
            print('Warming: the aggregate mode will be ignore for DSA data!')
        data_df = data_df[foldchange_threshold < abs(data_df['log2FoldChange'])]
        data_df = data_df[abs(data_df['log2FoldChange']) < float('inf')]
        data_df = data_df[data_df.iloc[:, 12:].apply(lambda x: np.any(x > avg_threshold) ,axis=1)]
        data_df2 = data_df['ReadCount1'].str.split(',', expand=True).astype(float)
        data_df2.columns = samples
        data_df2.index = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
        data_df = data_df2.apply(lambda x: np.log2(x + 1e-2))
    return data_df


def generate_heatmaps(data_df, original_columns, conditions, sample_cond_dict,
                      method, metric, prefix, aggregate, unsupervised, out_dir):
    num = data_df.shape[0]
    if num > 1000:
        print('Warming: number of rows > 1000, however, Jutils only allow 1000 rows for heatmap!')
        data_df = data_df.iloc[:1000, :]

    mask = data_df.isnull()
    data_df = data_df.fillna(0)

    clustermapParams = {
        'square':False # Tried to set this to True before. Don't: the dendograms do not scale well with it.
    }
    figureWidth = data_df.shape[1] / 4
    figureHeight= data_df.shape[0] / 4
    if figureHeight < 15:
        figureHeight = 15
    if figureWidth < 15:
        figureWidth = 15

    strings = [prefix]
    if unsupervised:
        strings.append('unsupervised')

    if aggregate:
        strings.append('aggregated')

    legend_title = 'PSI' if 'dPSI' in original_columns else 'log2(exp)'

    figure = sns.clustermap(data_df, cmap="RdBu_r", col_cluster=False, z_score=0, vmin=-5, vmax=5,
                            metric=metric, method=method, mask=mask,
                            yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight), **clustermapParams)

    figure.ax_heatmap.set_facecolor("lightyellow")
    set_xtick_text_colors(figure, sample_cond_dict, conditions)
    figure.ax_cbar.set_title(f'{legend_title}_Z', fontsize=16)
    figure.savefig(out_dir / '_'.join(filter(None, strings + ['clustermap_Z.png'])))
    plt.close()

    figure = sns.clustermap(data_df, cmap="RdBu_r", z_score=0, vmin=-5, vmax=5,
                            metric=metric, method=method, mask=mask,
                            yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight), **clustermapParams)

    figure.ax_heatmap.set_facecolor("lightyellow")
    set_xtick_text_colors(figure, sample_cond_dict, conditions)
    figure.ax_cbar.set_title(f'{legend_title}_Z', fontsize=16)
    figure.savefig(out_dir / '_'.join(filter(None, strings + ['clustermap2_Z.png'])))
    plt.close()

    if 'dPSI' in original_columns:
        figure = sns.clustermap(data_df, cmap=sns.cm.rocket_r, col_cluster=False,
                                metric=metric, method=method, mask=mask,
                                yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight), **clustermapParams)

        figure.ax_heatmap.set_facecolor("lightyellow")
        set_xtick_text_colors(figure, sample_cond_dict, conditions)
        figure.ax_cbar.set_title(legend_title, fontsize=16)
        figure.savefig(out_dir / '_'.join(filter(None, strings + ['clustermap.png'])))
        plt.close()

        figure = sns.clustermap(data_df, cmap=sns.cm.rocket_r,
                                metric=metric, method=method, mask=mask,
                                yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight), **clustermapParams)

        figure.ax_heatmap.set_facecolor("lightyellow")
        set_xtick_text_colors(figure, sample_cond_dict, conditions)
        figure.ax_cbar.set_title(legend_title, fontsize=16)
        figure.savefig(out_dir / '_'.join(filter(None, strings + ['clustermap2.png'])))
        plt.close()


def get_preset_palette():
    # palette = "#ff0000", "#00ff00", "#0000ff", "#000000"
    palette = ['#e6194B', '#000075', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#3cb44b', '#a9a9a9', '#000000']
    return palette


def set_xtick_text_colors(figure, sample_cond_dict, conditions):
    palette = get_preset_palette()
    n = len(palette)
    cond_color_dict = {cond: palette[i % n] for i, cond in enumerate(set(conditions))}
    for tick_label in figure.ax_heatmap.axes.get_xticklabels():
        cond = sample_cond_dict[tick_label.get_text()]
        tick_label.set_color(cond_color_dict[cond])

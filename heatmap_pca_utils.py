import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def check_thresholds(p_value_threshold, q_value_threshold, dpsi_threshold):
    if not(0 <= p_value_threshold <= 1):
        raise Exception('p-value threshold must in range [0, 1]!')
    if not(0 <= q_value_threshold <= 1):
        raise Exception('q-value threshold must in range [0, 1]!')
    if not(0 <= dpsi_threshold <= 1):
        raise Exception('dpsi threshold must in range [0, 1]!')


def print_warning(original_columns, p_value, q_value, dpsi, foldchange, avg, unsupervised, aggregate, top):
    if unsupervised:
        if p_value != 0.05:
            print('Warning: unsupervised mode, p-value option will be ignored!')
        if q_value != 1.:
            print('Warning: unsupervised mode, q-value option will be ignored!')
        if dpsi != 0.05:
            print('Warning: unsupervised mode, dpsi option will be ignored!')
        if foldchange != 0:
            print('Warning: unsupervised mode, fold-change option will be ignored!')
        if 'dPSI' in original_columns:
            if avg != 0:
                print('Warning: processing DSR data, avg option will be ignored!')

    else:
        if 'dPSI' in original_columns:
            if avg != 0:
                print('Warning: processing DSR data, avg option will be ignored!')
            if foldchange != 0:
                print('Warning: processing DSR data, fold-change option will be ignored!')
        elif 'log2FoldChange' in original_columns:
            if dpsi != 0.05:
                print('Warning: processing DSA data, dpsi option will be ignored!')
        if top != 100:
            print('Warning: top option only works in unsupervised mode and will be ignored!')

def parse_gene_list(gene_list_file):
    gene_list=[]
    with open(gene_list_file) as fi:
        for line in fi:
            if ' ' in line:
                sys.exit('Invalid format: space. Exiting...')
            else:
                gene_list.append(line.strip())
    return gene_list

def plot_heatmap_pca(file, meta_file, out_dir, p_value_threshold, q_value_threshold,
                 dpsi_threshold, foldchange_threshold, avg_threshold,
                 unsupervised, aggregate, prefix, top, pdf, gene_list_file, plot_type, method=None, metric=None, color_shape_col=None, label_point=None):
    check_thresholds(p_value_threshold, q_value_threshold, dpsi_threshold)

    samples = []
    conditions = []
    sample_cond_dict = {}
    with open(meta_file, 'r') as f:
        for line in f:
            sample, cond = line.strip().split('\t')[:2]
            samples.append(sample)
            conditions.append(cond)
            sample_cond_dict[sample] = cond

    data_df = pd.read_csv(file, sep='\t', comment='#', dtype={'GeneName': str})
    if gene_list_file:
        gene_list=parse_gene_list(gene_list_file)
        data_df=data_df[data_df['GeneName'].isin(gene_list)]
    data_df['GeneName'] = data_df['GeneName'].apply(truncate_gene_name)

    original_columns = data_df.columns
    if np.all(data_df['ReadCount1'] == '.') and np.all(data_df['PSI'] == '.'):
        raise Exception("Column 'ReadCount1' and 'PSI' don't contain values! This could happen if the TSV file is generated from MAJIQ outputs.")
        sys.exit('-1')

    print_warning(original_columns, p_value_threshold, q_value_threshold, dpsi_threshold, foldchange_threshold,
                  avg_threshold, unsupervised, aggregate, top)

    if plot_type=='pca' and 'log2FoldChange' in original_columns:
        for col in ['ReadCount1'] + data_df.columns[-len(set(conditions)):].to_list():
            data_df[col]=data_df[col].str.replace('None','0')

    if unsupervised:
        data_df = process_data_unsupervised(data_df, samples, original_columns, avg_threshold, aggregate, top)
    else:
        data_df = process_data_supervised(data_df, samples, sample_cond_dict, conditions, original_columns,
                                          p_value_threshold, q_value_threshold, dpsi_threshold, foldchange_threshold,
                                          avg_threshold, aggregate, plot_type)
    if plot_type=='heatmap':
        generate_heatmaps(data_df, original_columns, conditions, sample_cond_dict,
                          method, metric, prefix, aggregate, unsupervised, out_dir, pdf)
    elif plot_type=='pca':
        generate_pca(data_df, original_columns, meta_file, prefix, aggregate, unsupervised, out_dir, pdf, color_shape_col, label_point)


def truncate_gene_name(string):
    gene_names = string.split(',')
    if len(gene_names) > 3:
        return ','.join(gene_names[:3])
    return string


def process_group(x):
    start, end = float('inf'), 0
    for label in x['FeatureLabel']:
        if ',' in label or label.count(':') == 2:
            # rmats
            # chrX:111727445,111727613-111727770,111728184
            # chr20:58909349-58909423:58909520-5890957
            if len(x['FeatureLabel']) <= 1:
                return f"{x['GeneName'].iloc[0]}_{label}"
            else:
                l = len(x['FeatureLabel'])
                raise ValueError(f'rMATs data has {l} records per group')
        # non-rmats
        _chr, coord = label.split(':')
        _start, _end = (int(v) for v in coord.split('-'))
        start, end = min(start, _start), max(end, _end)
    return f"{x['GeneName'].iloc[0]}_{_chr}:{start}-{end}"

def process_data_unsupervised(data_df, samples, original_columns, avg_threshold, aggregate, top_num):
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
            data_df2 = data_df['PSI'].str.split(',', expand=True).replace('NA', '0').astype(float)
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
            data_df2 = data_df['PSI'].str.split(',', expand=True).replace('NA', '0').astype(float)
            data_df2.columns = samples

            # data_df2['variance'] = data_df2.mean(axis=1) / data_df2.var(axis=1)
            data_df2['variance'] = data_df2.sub(data_df2.mean(axis=1), axis=0).abs().mean(axis=1)
            if (data_df['FeatureType'] == 'intron').all():
                data_df2['index'] = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
            else:
                data_df2['index'] = data_df[['GeneName', 'FeatureLabel', 'FeatureType']].agg('_'.join, axis=1)

            groups = data_df2.groupby(['index'])
            data_df = groups.apply(lambda x: x.iloc[x['variance'].argmax()])
            data_df = data_df.sort_values(by=['variance'], ascending=False)
            data_df = data_df.drop(columns=['variance', 'index'])

    elif 'log2FoldChange' in original_columns:
        if aggregate:
            print('Warning: the aggregate mode will be ignore for DSA data!')
        data_df = data_df[data_df.iloc[:, 12:].apply(lambda x: np.any(x >= avg_threshold) ,axis=1)]
        data_df2 = data_df['ReadCount1'].str.split(',', expand=True).replace('NA', '0').astype(float)
        data_df2.columns = samples
        data_df2['variance'] = data_df2.var(axis=1) / data_df2.mean(axis=1)
        data_df2 = data_df2.apply(lambda x: np.log2(x + 1e-2))
        # data_df2['variance'] = data_df2.mean(axis=1) / data_df2.var(axis=1)
        if (data_df['FeatureType'] == 'intron').all():
            data_df2.index = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
        else:
            data_df2.index = data_df[['GeneName', 'FeatureLabel', 'FeatureType']].agg('_'.join, axis=1)
        data_df = data_df2
        data_df = data_df.sort_values(by=['variance'], ascending=False)
        data_df = data_df.drop(columns=['variance'])

    num = data_df.shape[0]
    if num > top_num:
        print(f'There are {num} events and the top {top_num} events shown (use "--top" option to change the number of events)!')
        num = min(top_num, data_df.shape[0])
    data_df = data_df.iloc[:num, :]
    return data_df


def process_data_supervised(data_df, samples, sample_cond_dict, conditions, original_columns,
                            p_value_threshold, q_value_threshold, dpsi_threshold,
                            foldchange_threshold, avg_threshold, aggregate, plot_type):
    data_df = data_df[data_df['p-value'] <= p_value_threshold]
    data_df = data_df[data_df['q-value'] <= q_value_threshold]
    if plot_type=='heatmap':
        annotated_gene_only=True
        max_dpsi_intron_only=True
    elif plot_type=='pca':
        annotated_gene_only=False
        max_dpsi_intron_only=False
    if annotated_gene_only:
        if np.any(data_df['GeneName'] != '.'):
            data_df = data_df[data_df['GeneName'] != '.']


    if 'dPSI' in original_columns:
        if data_df.shape[1] > 14 and aggregate:
            print("Multi-way comparison doesn't support aggregate mode, aggregate mode disabled!")
            aggregate = False

        if aggregate:
            selected_groups = data_df[abs(data_df['dPSI']) >= dpsi_threshold]['GroupID'].drop_duplicates()
            data_df = data_df[data_df['GroupID'].isin(selected_groups)]

            data_df2 = data_df['PSI'].str.split(',', expand=True).replace('NA', '0').astype(float)
            if data_df2.empty:
                sys.exit('Result set is empty. Exiting...')
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
            data_df = data_df[abs(data_df['dPSI']) >= dpsi_threshold]
            data_df2 = data_df['PSI'].str.split(',', expand=True).replace('NA', '0').astype(float)
            if data_df2.empty:
                sys.exit('Result set is empty. Exiting...')
            data_df2.columns = samples
            if max_dpsi_intron_only==False:
                data_df=data_df2
            else:
                if (data_df['FeatureType'] == 'intron').all():
                    data_df2['index'] = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
                else:
                    data_df2['index'] = data_df[['GeneName', 'FeatureLabel', 'FeatureType']].agg('_'.join, axis=1)
                data_df2['dPSI'] = np.abs(data_df['dPSI'])

                groups = data_df2.groupby(['index'])
                data_df = groups.apply(lambda x: x.iloc[x['dPSI'].argmax()])
                data_df = data_df.sort_values(by=['dPSI'], ascending=False)
                data_df = data_df.drop(columns=['dPSI', 'index'])

    elif 'log2FoldChange' in original_columns:
        if aggregate:
            print('Warning: the aggregate mode will be ignore for DSA data!')
        data_df = data_df[foldchange_threshold <= abs(data_df['log2FoldChange'])]
        data_df = data_df[abs(data_df['log2FoldChange']) <= float('inf')]
        data_df = data_df[data_df.iloc[:, 12:].astype(float).apply(lambda x: np.any(x >= avg_threshold) ,axis=1)]
        data_df2 = data_df['ReadCount1'].str.split(',', expand=True).replace('NA', '0').astype(float)
        data_df2.columns = samples
        if (data_df['FeatureType'] == 'intron').all():
            data_df2.index = data_df[['GeneName', 'FeatureLabel']].agg('_'.join, axis=1)
        else:
            data_df2.index = data_df[['GeneName', 'FeatureLabel', 'FeatureType']].agg('_'.join, axis=1)
        data_df = data_df2.apply(lambda x: np.log2(x + 1e-2))
    return data_df


def generate_heatmaps(data_df, original_columns, conditions, sample_cond_dict,
                      method, metric, prefix, aggregate, unsupervised, out_dir, pdf):
    num = data_df.shape[0]
    if num > 1000:
        print('Warning: number of rows > 1000, however, Jutils only allow 1000 rows for heatmap!')
        data_df = data_df.iloc[:1000, :]
    elif num < 2:
        print('Warning: number of rows < 2, data are not enough to generate heatmaps. Skipping...')
        return

    mask = data_df.isnull()
    data_df = data_df.fillna(0)

    clustermapParams = {
        'square':False # Tried to set this to True before. Don't: the dendograms do not scale well with it.
    }
    figureWidth = data_df.shape[1] / 4
    figureHeight= data_df.shape[0] / 4
    num_row=data_df.shape[0]

    # step function for dendrogram_height_ratio
    step_func_upper_bound=0.1
    step_func_lower_bound=0.03
    if num_row<=100:
        dendrogram_height_ratio=step_func_upper_bound
    elif num_row>100 and num_row<=500:
        dendrogram_height_ratio=step_func_upper_bound-(step_func_upper_bound-step_func_lower_bound)*(num_row-100)/(500-100)
    elif num_row>500:
        dendrogram_height_ratio=step_func_lower_bound
    #dendrogram_height_ratio=0.0008*data_df.shape[1]
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

    format = 'pdf' if pdf else 'png'
    figure = sns.clustermap(data_df, cmap="RdBu_r", col_cluster=False, z_score=0, vmin=-5, vmax=5,
                            metric=metric, method=method, mask=mask,
                            yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight),dendrogram_ratio=(0.2,dendrogram_height_ratio),cbar_pos=(0.02, 1-dendrogram_height_ratio, .02, .03), **clustermapParams)
    
    #  put color bar at the right
    #figure.fig.subplots_adjust(right=0.7)
    #figure.ax_cbar.set_position((0.02, .95, .02, .04))
    
    figure.ax_heatmap.set_facecolor("lightyellow")
    set_xtick_text_colors(figure, sample_cond_dict, conditions)
    figure.ax_cbar.set_title(f'{legend_title}_Z', fontsize=16)
    figure.savefig(out_dir / '_'.join(filter(None, strings + [f'clustermap_Z.{format}'])))
    plt.close()

    figure = sns.clustermap(data_df, cmap="RdBu_r", z_score=0, vmin=-5, vmax=5,
                            metric=metric, method=method, mask=mask,
                            yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight),dendrogram_ratio=(0.2,dendrogram_height_ratio),cbar_pos=(0.02, 1-dendrogram_height_ratio, .02, .03), **clustermapParams)

    figure.ax_heatmap.set_facecolor("lightyellow")
    set_xtick_text_colors(figure, sample_cond_dict, conditions)
    figure.ax_cbar.set_title(f'{legend_title}_Z', fontsize=16)
    figure.savefig(out_dir / '_'.join(filter(None, strings + [f'clustermap2_Z.{format}'])))
    plt.close()

    if 'dPSI' in original_columns:
        figure = sns.clustermap(data_df, cmap=sns.cm.rocket_r, col_cluster=False,
                                metric=metric, method=method, mask=mask,
                                yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight),dendrogram_ratio=(0.2,dendrogram_height_ratio),cbar_pos=(0.02, 1-dendrogram_height_ratio, .02, .03), **clustermapParams)

        figure.ax_heatmap.set_facecolor("lightyellow")
        set_xtick_text_colors(figure, sample_cond_dict, conditions)
        figure.ax_cbar.set_title(legend_title, fontsize=16)
        figure.savefig(out_dir / '_'.join(filter(None, strings + [f'clustermap.{format}'])))
        plt.close()

        figure = sns.clustermap(data_df, cmap=sns.cm.rocket_r,
                                metric=metric, method=method, mask=mask,
                                yticklabels=1, xticklabels=1, figsize=(figureWidth, figureHeight),dendrogram_ratio=(0.2,dendrogram_height_ratio),cbar_pos=(0.02, 1-dendrogram_height_ratio, .02, .03), **clustermapParams)

        figure.ax_heatmap.set_facecolor("lightyellow")
        set_xtick_text_colors(figure, sample_cond_dict, conditions)
        figure.ax_cbar.set_title(legend_title, fontsize=16)
        figure.savefig(out_dir / '_'.join(filter(None, strings + [f'clustermap2.{format}'])))
        plt.close()

def generate_pca(data_df, original_columns, meta_file, prefix, aggregate, unsupervised, out_dir, pdf, color_shape_col, label_point):

    from sklearn.decomposition import PCA
    import matplotlib

    def categorical_cmap(nc, nsc, cmap="tab10", continuous=False):
        if nc > plt.get_cmap(cmap).N:
            raise ValueError("Too many categories for colormap.")
        if continuous:
            ccolors = plt.get_cmap(cmap)(np.linspace(0,1,nc))
        else:
            ccolors = plt.get_cmap(cmap)(np.arange(nc, dtype=int))
        cols = np.zeros((nc*nsc, 3))
        for i, c in enumerate(ccolors):
            chsv = matplotlib.colors.rgb_to_hsv(c[:3])
            arhsv = np.tile(chsv,nsc).reshape(nsc,3)
            arhsv[:,1] = np.linspace(chsv[1],0.25,nsc)
            arhsv[:,2] = np.linspace(chsv[2],1,nsc)
            rgb = matplotlib.colors.hsv_to_rgb(arhsv)
            cols[i*nsc:(i+1)*nsc,:] = rgb
        cmap = matplotlib.colors.ListedColormap(cols)
        return cmap

    condition_col_index,group_col_index=color_shape_col.split(',')
    meta_df=pd.read_csv(meta_file,sep='\t',header=None)
    #  sort condition col
    condition_col_index=int(condition_col_index)-2
    meta_df=meta_df.sort_values([condition_col_index+1])
    label_df=meta_df.iloc[:,1:]
    labels_unique = label_df.drop_duplicates()
    targets = labels_unique.values


    conditions=labels_unique.iloc[:,condition_col_index]
    conditions_unique=conditions.unique()
    meta_col_n=meta_df.shape[1]
    if meta_col_n>2:
        group_col_index=int(group_col_index)-2
        groups_unique=labels_unique.iloc[:,group_col_index].unique()
        n_groups=len(groups_unique)
    else:
        n_groups=1
    n_conditions=len(conditions_unique)
    n_targets = len(targets)

    anno_text_len=20

    if 'dPSI' in original_columns:
        data_df=data_df.dropna()
        title='DSR'
    elif 'log2FoldChange' in original_columns:
        data_df=data_df.fillna(0)
        title='DSA'

    pca = PCA(n_components=3)
    principalComponents = pca.fit_transform(data_df.T)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['PC 1', 'PC 2', 'PC 3'])
    finalDf = pd.concat([principalDf, meta_df], axis = 1)
    for x,y in [(1,2),(1,3),(2,3)]:
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(f'PC {x}', fontsize = 15)
        ax.set_ylabel(f'PC {y}', fontsize = 15)
        ax.set_title(f'{title} 2-component PCA', fontsize = 20)

        colors=categorical_cmap(n_conditions, conditions.value_counts().max(), cmap="tab10").colors

        if meta_col_n>2:
            markers=dict()
            for i in range(1, n_groups+1):
                markers[groups_unique[i-1]]=(2+i, 1+i%2, i/n_targets*90.0)

        for i in range(n_targets):
            indicesToKeep = (finalDf.iloc[:,4:] == targets[i]).sum(axis=1) == meta_col_n-1
            x_pos=finalDf.loc[indicesToKeep, f'PC {x}']
            y_pos=finalDf.loc[indicesToKeep, f'PC {y}']
            if meta_col_n>2:
                ax.scatter(x_pos, y_pos, color=colors[i], marker=markers[targets[i][group_col_index]],s = 100)
            else:
                ax.scatter(x_pos, y_pos, color=colors[i],s = 50)
            if label_point==True:
                for j,row in finalDf.loc[indicesToKeep, [f'PC {x}',f'PC {y}',0]].iterrows():
                    ax.annotate(row[0][:anno_text_len],(row[f'PC {x}'], row[f'PC {y}']))
        legends=[]
        for target in targets:
            legends.append(' '.join([str(_) for _ in target]))

        ax.legend(legends)
        ax.grid()

        strings = [prefix]
        if unsupervised:
            strings.append('unsupervised')
        if aggregate:
            strings.append('aggregated')

        format = 'pdf' if pdf else 'png'
        fig.savefig(out_dir/'_'.join(filter(None, strings + [f'pca.pc{x}-{y}.{format}'])))


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

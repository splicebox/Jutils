from collections import defaultdict, OrderedDict
import numpy as np
import gzip, os

from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def convert_majiq_results(data_dir, out_dir):
    infile = None
    for file in os.listdir(data_dir):
        if file.endswith(".tsv"):
            infile = file
            break

    if not infile:
        raise FileNotFoundError(f'Error: the .tsv file not found in {data_dir}!')

    with open(data_dir / infile, 'r') as f:
        lines = f.readlines()

    intron_info_dict = OrderedDict()
    conds = [f'psi({c[:-7]})' for c in lines[0].strip().split('\t')[6:8]]
    for line in lines[1:]:
        items = line.strip().split('\t')
        gene_name = items[0]
        lsv_id = items[2]
        _chr, strand = items[15], items[16]

        dpsis = items[3].split(';')
        intron_coords = [v.split('-') for v in items[17].split(';')]
        cond1_psis = items[6].split(';')
        cond2_psis = items[7].split(';')

        i = 1
        for dpsi, (start, end), psi1, psi2 in zip(dpsis, intron_coords, cond1_psis, cond2_psis):
            start, end = int(start), int(end)
            dpsi = float(dpsi)
            intron_info_dict[(_chr, start, end, lsv_id)] = f'{gene_name}\t{lsv_id}\ti{i:03d}\tintron\t{_chr}:{start}-{end}\t{strand}\t.\t.\t{dpsi:.6g}\t.\t.\t.\t{psi1}\t{psi2}'
            i += 1

    out_buffer = 'GeneName\tGroupID\tFeatureID\tFeatureType\tFeatureLabel\tstrand\tp-value\tq-value\tdPSI\tReadCount1\tReadCount2\tPSI\t' + '\t'.join(conds)
    file = out_dir / 'majiq_results.tsv'
    with open(file, 'w') as f:
        f.write('# majiq\n')
        f.write(out_buffer + '\n')
        for info in intron_info_dict.values():
            f.write(info + '\n')


def convert_leafcutter_results(data_dir, out_dir):
    def custom_divide(psi_str):
        a, b = [float(v) for v in psi_str.split('/')]
        if b == 0:
            return 'nan'
        return f'{a/b:.6g}'

    data_dir, out_dir = Path(data_dir), Path(out_dir)
    file = data_dir / 'leafcutter_ds_cluster_significance.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    cluster_info_dict = {}
    for line in lines[1:]:
        cluster, status, _, _, p_value, q_value, gene_names_str = line.strip().split('\t')
        if status == 'Success':
            _chr, cluster_id = cluster.split(':')
            gene_names_str = '.' if gene_names_str == 'NA' else gene_names_str
            p_value = f'{float(p_value):.6g}'
            q_value = f'{float(q_value):.6g}'
            cluster_info_dict[cluster_id] = (p_value, q_value, gene_names_str)

    file = data_dir / 'leafcutter_ds_effect_sizes.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    intron_info_dict = OrderedDict()
    intron_info_psis_dict = {}
    prev_cluster, i = None, 0
    conds = [f'psi({c})' for c in lines[0].strip().split('\t')[2:4]]
    for line in lines[1:]:
        # intron  logef   case    control deltapsi
        intron_info, _, psi1, psi2, dpsi = line.strip().split('\t')
        _chr, start, end, cluster_id = intron_info.split(':')
        strand = cluster_id.split('_')[-1]
        strand = strand if strand != 'NA' else '.'
        p_value, q_value, gene_names_str = cluster_info_dict[cluster_id]
        dpsi = f'{float(dpsi):.6g}'
        i = 1 if prev_cluster != cluster_id else i + 1
        prev_cluster = cluster_id
        intron_info_dict[intron_info] = f'{gene_names_str}\t{cluster_id}\ti{i:03d}\tintron\t{_chr}:{start}-{end}\t{strand}\t{p_value}\t{q_value}\t{dpsi}'
        intron_info_psis_dict[intron_info] = (psi1, psi2)

    file = data_dir / 'results_perind_numers.counts.gz'
    with gzip.open(file, 'rb') as f:
        lines = f.readlines()

    for line in lines:
        items = line.decode('utf-8').strip().split(' ')
        intron_info = items[0]
        if intron_info in intron_info_dict:
            intron_info_dict[intron_info] += f"\t{','.join(items[1:])}\t."

    file = data_dir / 'results_perind.counts.gz'
    with gzip.open(file, 'rb') as f:
        lines = f.readlines()

    for line in lines:
        items = line.decode('utf-8').split(' ')
        intron_info = items[0]
        if intron_info in intron_info_dict:
            psis = [custom_divide(v) for v in items[1:]]
            intron_info_dict[intron_info] += f"\t{','.join(psis)}\t" + '\t'.join(intron_info_psis_dict[intron_info])

    out_buffer = 'GeneName\tGroupID\tFeatureID\tFeatureType\tFeatureLabel\tstrand\tp-value\tq-value\tdPSI\tReadCount1\tReadCount2\tPSI\t' + '\t'.join(conds)
    file = out_dir / 'leafcutter_results.tsv'
    with open(file, 'w') as f:
        f.write('# leafcutter\n')
        f.write(out_buffer + '\n')
        for info in intron_info_dict.values():
            f.write(info + '\n')


def convert_rmats_results(data_dir, out_dir):
    data_dir, out_dir = Path(data_dir), Path(out_dir)
    for file_type in ['ReadsOnTargetAndJunctionCounts', 'JunctionCountOnly']:
        out_buffer = 'GeneName\tGroupID\tFeatureID\tFeatureType\tFeatureLabel\tstrand\tp-value\tq-value\tdPSI\tReadCount1\tReadCount2\tPSI\tpsi(cond1)\tpsi(cond2)\n'
        #ID    GeneID   geneSymbol  chr  strand  exonStart_0base  exonEnd   upstreamES  upstreamEE  downstreamES  downstreamEE  ID    IC_SAMPLE_1
        file = data_dir / f'SE.MATS.{file_type}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()

        for line in lines[1:]:
            fid, _, gene_name, _chr, strand, es, ee, ues, uee, des, dee, _, ic1, sc1, ic2, sc2, _, _, p_value, q_value, icl1, icl2, dpsi = line.strip().split('\t')
            label = f'{_chr}:{uee},{es}-{ee},{des}'
            p_value, q_value, dpsi = f'{float(p_value):.6g}', f'{float(q_value):.6g}', f'{float(dpsi):.6g}'
            gene_name = gene_name[1:-1]
            out_buffer += f'{gene_name}\t.\t{fid}\tSE\t{label}\t{strand}\t{p_value}\t{q_value}\t{dpsi}\t{ic1},{ic2}\t{sc1},{sc2}\t{icl1},{icl2}\t.\t.\n'

        file = data_dir / f'RI.MATS.{file_type}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()

        for line in lines[1:]:
            fid, _, gene_name, _chr, strand, _, _, ues, uee, des, dee, _, ic1, sc1, ic2, sc2, _, _, p_value, q_value, icl1, icl2, dpsi = line.strip().split('\t')
            label = f'{_chr}:{ues}-{uee}:{des}-{dee}'
            p_value, q_value, dpsi = f'{float(p_value):.6g}', f'{float(q_value):.6g}', f'{float(dpsi):.6g}'
            gene_name = gene_name[1:-1]
            out_buffer += f'{gene_name}\t.\t{fid}\tRI\t{label}\t{strand}\t{p_value}\t{p_value}\t{dpsi}\t{ic1},{ic2}\t{sc1},{sc2}\t{icl1},{icl2}\t.\t.\n'

        file = data_dir / f'MXE.MATS.{file_type}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()

        for line in lines[1:]:
            fid, _, gene_name, _chr, strand, es1, ee1, es2, ee2, ues, uee, des, dee, _, ic1, sc1, ic2, sc2, _, _, p_value, q_value, icl1, icl2, dpsi = line.strip().split('\t')
            label = f'{_chr}:{uee},{es1}-{ee1}:{es2}-{ee2},{des}'
            p_value, q_value, dpsi = f'{float(p_value):.6g}', f'{float(q_value):.6g}', f'{float(dpsi):.6g}'
            gene_name = gene_name[1:-1]
            out_buffer += f'{gene_name}\t.\t{fid}\tMXE\t{label}\t{strand}\t{p_value}\t{p_value}\t{dpsi}\t{ic1},{ic2}\t{sc1},{sc2}\t{icl1},{icl2}\t.\t.\n'

        file = data_dir / f'A5SS.MATS.{file_type}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()

        for line in lines[1:]:
            fid, _, gene_name, _chr, strand, les, lee, ses, see, fes, fee, _, ic1, sc1, ic2, sc2, _, _, p_value, q_value, icl1, icl2, dpsi = line.strip().split('\t')
            label = f'{_chr}:{les}-{lee}:{ses}-{see},{fes}'
            p_value, q_value, dpsi = f'{float(p_value):.6g}', f'{float(q_value):.6g}', f'{float(dpsi):.6g}'
            gene_name = gene_name[1:-1]
            out_buffer += f'{gene_name}\t.\t{fid}\tA5SS\t{label}\t{strand}\t{p_value}\t{p_value}\t{dpsi}\t{ic1},{ic2}\t{sc1},{sc2}\t{icl1},{icl2}\t.\t.\n'

        file = data_dir / f'A3SS.MATS.{file_type}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()

        for line in lines[1:]:
            fid, _, gene_name, _chr, strand, les, lee, ses, see, fes, fee, _, ic1, sc1, ic2, sc2, _, _, p_value, q_value, icl1, icl2, dpsi = line.strip().split('\t')
            label = f'{_chr}:{fee},{les}-{lee}:{ses}-{see}'
            p_value, q_value, dpsi = f'{float(p_value):.6g}', f'{float(q_value):.6g}', f'{float(dpsi):.6g}'
            gene_name = gene_name[1:-1]
            out_buffer += f'{gene_name}\t.\t{fid}\tA3SS\t{label}\t{strand}\t{p_value}\t{p_value}\t{dpsi}\t{ic1},{ic2}\t{sc1},{sc2}\t{icl1},{icl2}\t.\t.\n'

        file = out_dir / f'rmats_{file_type}_results.tsv'
        with open(file, 'w') as f:
            f.write('# rmats\n')
            f.write(out_buffer)


def convert_mntjulip_results(data_dir, out_dir):
    convert_mntjulip_DSR_results(data_dir, out_dir)
    convert_mntjulip_DSA_results(data_dir, out_dir)


def convert_mntjulip_DSR_results(data_dir, out_dir):
    def strs_to_floats(counts):
        return [float(c) for c in counts.split(',')]

    def custom_sum(introns, intron_counts_dict):
        array = [strs_to_floats(intron_counts_dict[i]) for i in introns]
        return np.sum(array, axis=0)

    def custom_divide(counts, sums):
        counts = strs_to_floats(counts)
        res = (counts / sums).tolist()
        res = ','.join([f'{i:.6g}' for i in res])
        return res

    def max_dpsi(items, indices):
        max_dpsi = 0
        for k, i in enumerate(indices):
            for j in indices[k+1:]:
                psi1, psi2 = float(items[i]), float(items[j])
                dpsi = psi2 - psi1
                if len(indices) == 2:
                    return dpsi
                max_dpsi = max(max_dpsi, abs(dpsi))
        return max_dpsi

    data_dir, out_dir = Path(data_dir), Path(out_dir)
    file = data_dir / 'diff_spliced_groups.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    group_info_dict = {}
    group_id_name_dict = {}
    for line in lines[1:]:
        # group_id  chrom  loc  strand  gene_name  structure  llr  p_value   q_value
        group_id, _chr, loc, strand, gene_names_str, _, _, p_value, q_value = line.strip().split('\t')
        p_value, q_value = f'{float(p_value):.6g}', f'{float(q_value):.6g}'
        group_id_name_dict[group_id] = f"{gene_names_str}_{_chr}_{loc}"
        group_info_dict[group_id] = (p_value, q_value, gene_names_str)

    file = data_dir / 'diff_spliced_introns.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    intron_info_dict = OrderedDict()
    group_introns_dict = defaultdict(set)
    prev_group, i = None, 0
    intron_info_psis_dict = {}
    items = lines[0].strip().split('\t')
    indices = [i for i in range(len(items)) if items[i].startswith('psi')]
    conds = [items[i] for i in indices]
    for line in lines[1:]:
        # group_id  chrom  start  end  strand  gene_name  psi(Normal)  psi(Tumor)   delta_psi
        items = line.strip().split('\t')
        group_id, _chr, start, end, strand = items[0], items[1], items[2], items[3], items[4]
        p_value, q_value, gene_names_str = group_info_dict[group_id]
        intron = f'{_chr}:{start}-{end}'
        intron_info = f'{_chr}:{start}-{end}:{group_id}'
        i = 1 if prev_group != group_id else i + 1
        prev_group = group_id
        dpsi = max_dpsi(items, indices)
        intron_info_dict[intron_info] = f'{gene_names_str}\t{group_id}\ti{i:03d}\tintron\t{intron}\t{strand}\t{p_value}\t{q_value}\t{dpsi:.6g}'
        group_introns_dict[group_id].add(intron)
        intron_info_psis_dict[intron_info] = [items[i] for i in indices]

    file = data_dir / 'intron_data.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    intron_counts_dict = {}
    for line in lines[1:]:
        # chrom  start end  strand  gene_name  status  read_counts(Normal) read_counts(Tumor)
        items = line.strip().split('\t')
        _chr, start, end, strand = items[:4]
        intron = f'{_chr}:{start}-{end}'
        intron_counts_dict[intron] = ','.join(items[6:])

    for group_id, introns in group_introns_dict.items():
        sums = custom_sum(introns, intron_counts_dict)
        for intron in introns:
            intron_info = f'{intron}:{group_id}'
            intron_info_dict[intron_info] += f"\t{intron_counts_dict[intron]}\t.\t{custom_divide(intron_counts_dict[intron], sums)}\t"
            intron_info_dict[intron_info] += '\t'.join(intron_info_psis_dict[intron_info])

    out_buffer = 'GeneName\tGroupID\tFeatureID\tFeatureType\tFeatureLabel\tstrand\tp-value\tq-value\tdPSI\tReadCount1\tReadCount2\tPSI\t' + '\t'.join(conds)
    file = out_dir / 'mntjulip_DSR_results.tsv'
    with open(file, 'w') as f:
        f.write('# mntjulip DSR\n')
        f.write(out_buffer + '\n')
        for info in intron_info_dict.values():
            f.write(info + '\n')


def convert_mntjulip_DSA_results(data_dir, out_dir):
    def max_fold_change(items, indices):
        max_fold_change = 0
        for k, i in enumerate(indices):
            for j in indices[k+1:]:
                m1, m2 = float(items[i]), float(items[j])
                fold_change = 0 if m2 == 0 else float('inf')
                if m1 > 0 and m2 != 0:
                    fold_change = np.log2(m2 / m1)
                if len(indices) == 2:
                    return fold_change
                if max_fold_change < abs(fold_change) < float('inf'):
                    max_fold_change = abs(fold_change)
        return max_fold_change

    file = data_dir / 'diff_introns.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    intron_info_dict = {}
    intron_means_dict = {}
    items = lines[0].strip().split('\t')
    indices = [i for i in range(len(items)) if items[i].startswith('avg')]
    conds = [items[i] for i in indices]
    for line in lines[1:]:
        # chrom start end strand gene_name status  llr     p_value q_value avg_read_counts(case)   avg_read_counts(control)
        items = line.strip().split('\t')
        _chr, start, end, strand, gene_names_str, status = items[:6]
        if status == 'TEST':
            p_value, q_value = items[7: 9]
            intron = f'{_chr}:{start}-{end}'
            fold_change = max_fold_change(items, indices)
            intron_info_dict[(intron, strand)] = f'{gene_names_str}\t.\t.\tintron\t{intron}\t{strand}\t{p_value}\t{q_value}\t{fold_change:.6g}\t'
            intron_means_dict[(intron, strand)] = [items[i] for i in indices]
    file = data_dir / 'intron_data.txt'
    with open(file, 'r') as f:
        lines = f.readlines()

    intron_counts_dict = {}
    for line in lines[1:]:
        # chrom  start end  strand  gene_name  status  read_counts(Normal) read_counts(Tumor)
        items = line.strip().split('\t')
        _chr, start, end, strand = items[:4]
        intron = f'{_chr}:{start}-{end}'
        if (intron, strand) in intron_info_dict:
            intron_info_dict[(intron, strand)] += ','.join(items[6:]) + '\t.\t.\t'
            intron_info_dict[(intron, strand)] += '\t'.join(intron_means_dict[(intron, strand)])

    out_buffer = f'GeneName\tGroupID\tFeatureID\tFeatureType\tFeatureLabel\tstrand\tp-value\tq-value\tlog2FoldChange\tReadCount1\tReadCount2\tPSI\t' + '\t'.join(conds)
    file = out_dir / 'mntjulip_DSA_results.tsv'
    with open(file, 'w') as f:
        f.write('# mntjulip DSA\n')
        f.write(out_buffer + '\n')
        for info in intron_info_dict.values():
            f.write(info + '\n')

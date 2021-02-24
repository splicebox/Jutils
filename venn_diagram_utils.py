import venn
import matplotlib.pyplot as plt
import argparse


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--p-value', type=float, default=0.05)
    parser.add_argument('--q-value', type=float, default=1.0)
    parser.add_argument('--dpsi', type=float, default=0.05)
    return parser


def plot_venn_diagram(list_file, out_dir, p_value_threhold, q_value_threhold, dpsi_threshold):
    files = []
    alias = []
    options = []
    with open(list_file, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            files.append(items[0])
            if len(items) > 1:
                alias.append(items[1])
            if len(items) > 2:
                options.append(items[2])

    if len(files) > 6:
        print('As most 6 files are allowed, the program will automatically choose the first 6 files!')
        files = files[:6]
        alias = alias[:6] if alias else []
        options = options[:6] if options else []

    if options:
        print('Options for each TSV file are provided, will overwrite global options!')
        parser = get_parser()

    names = []
    genes_list = []
    for i, file in enumerate(files):
        l_p_value_threhold, l_q_value_threhold, l_dpsi_threshold = p_value_threhold, q_value_threhold, dpsi_threshold
        if options and options[i]:
            arg = parser.parse_args(options[i][1:-1].split())
            l_p_value_threhold, l_q_value_threhold, l_dpsi_threshold = arg.p_value, arg.q_value, arg.dpsi
        genes = set()
        with open(file) as f:
            lines = f.readlines()

        names.append(lines[0].strip()[2:] if not alias else alias[i])
        for line in lines:
            if line.startswith('#') or line.startswith('GeneName'):
                continue

            items = line.split('\t')
            gene_names = items[0]
            p_value, q_value, dpsi = items[6: 9]
            p_value = 0 if p_value == '.' else float(p_value)
            q_value = 0 if q_value == '.' else float(q_value)
            dpsi = 1 if dpsi == '.' else float(dpsi)

            if p_value < l_p_value_threhold and q_value < l_q_value_threhold and abs(dpsi) > l_dpsi_threshold:
                for gene_name in gene_names.split('.'):
                    genes.add(gene_name)
        genes_list.append(genes)

    venn_method = getattr(venn, f'venn{len(files)}')
    labels = venn.get_labels(genes_list, fill=['number'])
    fig, ax = venn_method(labels, names=alias if alias else names)
    fig.patch.set_facecolor('white')
    fig.savefig(out_dir / 'venn_diagram.png')
    plt.close()

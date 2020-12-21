import venn
import matplotlib.pyplot as plt


def plot_venn_diagram(list_file, out_dir, p_value_threhold, q_value_threhold, dpsi_threshold):
    files = []
    alias = []
    with open(list_file, 'r') as f:
        for line in f:
            items = line.strip().split('\t')
            files.append(items[0])
            if len(items) > 1:
                alias.append(items[1])

    if len(files) > 6:
        print('As most 6 files are allowed, the program will automatically choose the first 6 files!')
        files = files[:6]
        alias = alias[:6] if alias else []

    names = []
    genes_list = []
    for i, file in enumerate(files):
        genes = set()
        with open(file) as f:
            lines = f.readlines()

        names.append(lines[0].strip()[2:] if not alias else alias[i])
        for line in lines[2:]:
            items = line.split('\t')
            gene_names = items[0]
            p_value, q_value, dpsi = items[6: 9]

            p_value = 0 if p_value == '.' else float(p_value)
            q_value = 0 if q_value == '.' else float(q_value)
            dpsi = 1 if dpsi == '.' else float(dpsi)
            if p_value < p_value_threhold and q_value < q_value_threhold and abs(dpsi) > dpsi_threshold:
                for gene_name in gene_names.split('.'):
                    genes.add(gene_name)
        genes_list.append(genes)

    venn_method = getattr(venn, f'venn{len(files)}')
    labels = venn.get_labels(genes_list, fill=['number'])
    fig, ax = venn_method(labels, names=alias if alias else names)
    fig.patch.set_facecolor('white')
    fig.savefig(out_dir / 'venn_diagram.png')
    plt.close()

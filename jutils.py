from pathlib import Path
import sys, argparse, textwrap

from convert_results_utils import convert_leafcutter_results, convert_rmats_results, convert_mntjulip_results, convert_majiq_results
from venn_diagram_utils import plot_venn_diagram
from heatmap_pca_utils import plot_heatmap_pca
from sashimi_utils import sashimi_plot_with_bams, sashimi_plot_without_bams


def get_arguments():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(help='commands', dest="command")
    parser_dict = {}

    r_parser = subparser.add_parser('convert-results')
    parser_dict['convert-results'] = r_parser
    r_parser.add_argument('--leafcutter-dir', type=str, help='the directory contains the leafcutter_ds_cluster_significance.txt, leafcutter_ds_effect_sizes.txt and *.counts.gz')
    r_parser.add_argument('--rmats-dir', type=str, help='the directory contains the *.ReadsOnTargetAndJunctionCounts.txt and *.JunctionCountOnly.txt.')
    r_parser.add_argument('--mntjulip-dir', type=str, help='the directory contains diff_spliced_introns.txt, diff_spliced_groups.txt, diff_introns.txt and intron_data.txt')
    r_parser.add_argument('--majiq-dir', type=str, help='the directory contains *.deltapsi.tsv')
    r_parser.add_argument('--out-dir', type=str, default='./out', help='specify the output directory')

    v_parser = subparser.add_parser('venn-diagram', help='')
    parser_dict['venn-diagram'] = v_parser
    v_parser.add_argument('--tsv-file-list', type=str, help='a file that contains path of tsv files')
    v_parser.add_argument('--p-value', type=float, default=0.05, help='cutoff for differential test p-value (default 1.0)')
    v_parser.add_argument('--q-value', type=float, default=1.0, help='cutoff for differential test q-value (default 1.0)')
    v_parser.add_argument('--out-dir', type=str, default='./out', help='specify the output directory')
    v_parser.add_argument('--dpsi', type=float, default=0.05, help='cutoff for delta PSI (Percent Splice In) (default 0.0)')
    v_parser.add_argument('--prefix', type=str, default='', help='add prefix to the output file')
    v_parser.add_argument('--pdf', action='store_true', default=False, help='generate figure(s) in .pdf format')

    h_parser = subparser.add_parser('heatmap', help='')
    parser_dict['heatmap'] = h_parser
    h_parser.add_argument('--tsv-file', type=str, help='a TSV file that contains the extracted results')
    h_parser.add_argument('--meta-file', type=str, help='a TAB separated file that contains the sample name and conditions')
    h_parser.add_argument('--p-value', type=float, default=0.05, help='cutoff for differential test p-value (default 0.05)')
    h_parser.add_argument('--q-value', type=float, default=1.0, help='cutoff for differential test q-value (default 1.0)')
    h_parser.add_argument('--dpsi', type=float, default=0.05, help='cutoff for delta PSI (Percent Splice In) (default 0.05)')
    h_parser.add_argument('--avg', type=float, default=0., help='cutoff for estimated read counts of DSA results (default 0.0)')
    h_parser.add_argument('--fold-change', type=float, default=0., help='cutoff for log2(fold-change) of DSA results (default 0.0)')
    h_parser.add_argument('--aggregate', action='store_true', default=False, help='show results at group level (one entry per group)')
    h_parser.add_argument('--unsupervised', action='store_true', default=False, help='show results in an unsupervised way')
    h_parser.add_argument('--out-dir', type=str, default='./out', help='specify the output directory')
    h_parser.add_argument('--prefix', type=str, default='', help='add prefix to the output file')
    h_parser.add_argument('--top', type=int, default=100, help='number of top most variable features to display, this option only works in the unsupervised mode (default 100)')
    h_parser.add_argument('--method', type=str, default='weighted',
                        help="clustering method. choose from 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' (default 'weighted')")
    h_parser.add_argument('--metric', type=str, default='cityblock',
                        help="distance metric for clustering. The distance function can be 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule' (default 'cityblock')")
    h_parser.add_argument('--pdf', action='store_true', default=False, help='generate figure(s) in .pdf format')
    h_parser.add_argument('--gene-list-file', type=str, default='', help='list of target genes (one gene per line without space) for heatmap')

    p_parser = subparser.add_parser('pca', help='')
    parser_dict['pca'] = p_parser
    p_parser.add_argument('--tsv-file', type=str, help='a TSV file that contains the extracted results')
    p_parser.add_argument('--meta-file', type=str, help='a TAB separated file that contains the sample name and conditions')
    p_parser.add_argument('--p-value', type=float, default=0.05, help='cutoff for differential test p-value (default 0.05)')
    p_parser.add_argument('--q-value', type=float, default=1.0, help='cutoff for differential test q-value (default 1.0)')
    p_parser.add_argument('--dpsi', type=float, default=0.05, help='cutoff for delta PSI (Percent Splice In) (default 0.05)')
    p_parser.add_argument('--avg', type=float, default=0., help='cutoff for estimated read counts of DSA results (default 0.0)')
    p_parser.add_argument('--fold-change', type=float, default=0., help='cutoff for log2(fold-change) of DSA results (default 0.0)')
    p_parser.add_argument('--aggregate', action='store_true', default=False, help='show results at group level (one entry per group)')
    p_parser.add_argument('--unsupervised', action='store_true', default=False, help='show results in an unsupervised way')
    p_parser.add_argument('--out-dir', type=str, default='./out', help='specify the output directory')
    p_parser.add_argument('--prefix', type=str, default='', help='add prefix to the output file')
    p_parser.add_argument('--top', type=int, default=100, help='number of top most variable features to display, this option only works in the unsupervised mode (default 100)')
    p_parser.add_argument('--pdf', action='store_true', default=False, help='generate figure(s) in .pdf format')
    p_parser.add_argument('--gene-list-file', type=str, default='', help='list of target genes (one gene per line without space) for heatmap')
    p_parser.add_argument('--color-shape-col', type=str, default='2,3', help='Colour and shape points by the 2 column indices in the meta file. The sample name column starts with index 1 (default: 2,3)')
    p_parser.add_argument('--label-point', action='store_true', default=False, help='label points with sample names from the meta file')

    s_parser = subparser.add_parser('sashimi', help='')
    parser_dict['sashimi'] = s_parser
    s_parser.add_argument('--bam-list', type=str, help='a BAM files list')
    s_parser.add_argument('--tsv-file', type=str, help='a TSV file that contains the extracted results')
    s_parser.add_argument('--meta-file', type=str, help='a TAB separated file that contains the sample name and conditions')
    s_parser.add_argument('--coordinate', type=str, help='select genomic range for visualization (e.g. chr1:123456-234567)')
    #  TODO gtf optional
    s_parser.add_argument('--gtf', type=str, help='a GTF file with gene annotations')
    s_parser.add_argument('--shrink', action='store_true', default=False, help='shrink long introns and exons')
    #  TODO filter min-coverage of average
    s_parser.add_argument("--min-coverage", type=int, default=0, help='minimum intron coverage')
    s_parser.add_argument("--group-id", type=str, help='Specify the group id')
    # s_parser.add_argument("--strand", default="NONE", type=str,
    #     help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
    s_parser.add_argument('--out-dir', type=str, default='./out', help='specify the output directory')
    s_parser.add_argument('--prefix', type=str, default='', help='add prefix to the output file')
    s_parser.add_argument('--pdf', action='store_true', default=False, help='generate figure(s) in .pdf format')


    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args(), parser_dict


def run_convert_results_module(args, parser_dict):
    
    if args.mntjulip_dir:
        convert_mntjulip_results(Path(args.mntjulip_dir), Path(args.out_dir))
    if args.leafcutter_dir:
        convert_leafcutter_results(Path(args.leafcutter_dir), Path(args.out_dir))
    if args.rmats_dir:
        convert_rmats_results(Path(args.rmats_dir), Path(args.out_dir))
    if args.majiq_dir:
        convert_majiq_results(Path(args.majiq_dir), Path(args.out_dir))
    if not (args.mntjulip_dir or args.leafcutter_dir or args.rmats_dir or args.majiq_dir):
        # raise Exception('Please specify at least a path to a program result folder!')
        parser_dict['convert-results'].print_help(sys.stderr)


def run_venn_diagram_module(args, parser_dict):
    if not args.tsv_file_list:
        # raise Exception('Please provide the list file that contains the path of the TSV result files!')
        parser_dict['venn-diagram'].print_help(sys.stderr)
    else:
        plot_venn_diagram(Path(args.tsv_file_list), Path(args.out_dir), args.p_value, args.q_value, args.dpsi, args.prefix, args.pdf)


def run_heatmap_module(args, parser_dict):
    if not args.tsv_file or not args.meta_file:
        # raise Exception('Please provide the list file that contains the path of the TSV result files!')
        parser_dict['heatmap'].print_help(sys.stderr)
    else:
        plot_heatmap_pca(Path(args.tsv_file), Path(args.meta_file), Path(args.out_dir), args.p_value,
                 args.q_value, args.dpsi, args.fold_change, args.avg, args.unsupervised,
                 args.aggregate, args.prefix, args.top, args.pdf, args.gene_list_file, plot_type='heatmap', method=args.method, metric=args.metric)

def run_pca_module(args, parser_dict):
    if not args.tsv_file or not args.meta_file:
        # raise Exception('Please provide the list file that contains the path of the TSV result files!')
        parser_dict['pca'].print_help(sys.stderr)
    else:
        plot_heatmap_pca(Path(args.tsv_file), Path(args.meta_file), Path(args.out_dir), args.p_value,
                 args.q_value, args.dpsi, args.fold_change, args.avg, args.unsupervised,
                 args.aggregate, args.prefix, args.top, args.pdf, args.gene_list_file, plot_type='pca', color_shape_col=args.color_shape_col, label_point=args.label_point)

def run_sashimi_module(args, parser_dict):
    if args.bam_list:
        sashimi_plot_with_bams(args.bam_list, args.coordinate, args.gtf, Path(args.out_dir),
                               args.prefix, args.shrink, 'NONE', args.min_coverage,
                               args.group_id, args.tsv_file, args.pdf)
    elif args.tsv_file and args.meta_file and args.gtf:
        sashimi_plot_without_bams(args.tsv_file, args.meta_file, args.gtf, args.group_id,
                                  Path(args.out_dir), args.prefix, args.shrink, args.min_coverage, args.pdf)
    else:
    #     raise Exception(textwrap.dedent('''\
    # Please provide one of the two set of files,
    # (1) A bam file list and a coordinate (e.g. chr1:12345-12356)
    # (2) A TSV file that are generated by the 'convert-results' module, a meta file of sample names and conditions, a GTF reference genome file
    # '''))
        parser_dict['sashimi'].print_help(sys.stderr)


def main():
    args, parser_dict = get_arguments()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.command == 'convert-results':
        run_convert_results_module(args, parser_dict)

    if args.command == 'venn-diagram':
        run_venn_diagram_module(args, parser_dict)

    if args.command == 'heatmap':
        run_heatmap_module(args, parser_dict)

    if args.command == 'pca':
        run_pca_module(args, parser_dict)

    if args.command == 'sashimi':
        run_sashimi_module(args, parser_dict)


if __name__ == "__main__":
    main()

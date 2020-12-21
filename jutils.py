from pathlib import Path
import sys, argparse, textwrap

from convert_results_utils import convert_leafcutter_results, convert_rmats_results, convert_mntjulip_results, convert_majiq_results
from venn_diagram_utils import plot_venn_diagram
from heatmap_utils import plot_heatmap
from sashimi_utils import sashimi_plot_with_bams, sashimi_polt_without_bams


def get_arguments():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(help='commands', dest="command")

    r_parser = subparser.add_parser('convert-results')
    r_parser.add_argument('--leafcutter-dir', type=str, help='The directory contains the leafcutter_ds_cluster_significance.txt, leafcutter_ds_effect_sizes.txt and *.counts.gz')
    r_parser.add_argument('--rmat-dir', type=str, help='The directory contains the *.ReadsOnTargetAndJunctionCounts.txt and *.JunctionCountOnly.txt.')
    r_parser.add_argument('--mntjulip-dir', type=str, help='The directory contains diff_spliced_introns.txt, diff_spliced_groups.txt, diff_introns.txt and intron_data.txt')
    r_parser.add_argument('--majiq-dir', type=str, help='The directory contains *.deltapsi.tsv')
    r_parser.add_argument('--out-dir', type=str, default='./out', help='The output directory')

    v_parser = subparser.add_parser('venn-diagram', help='')
    v_parser.add_argument('--tsv-file-list', type=str, help='A file that contains path of tsv files')
    v_parser.add_argument('--p-value', type=float, default=0.05, help='Provide a p-value cutoff (default 0.05)')
    v_parser.add_argument('--q-value', type=float, default=1.0, help='Provide a q-value cutoff (default 1.0)')
    v_parser.add_argument('--out-dir', type=str, default='./out', help='The output directory')
    v_parser.add_argument('--dpsi', type=float, default=0.05, help='Provide a |dPSI| cutoff (default 0.05)')
    v_parser.add_argument('--prefix', type=str, default='', help='Add prefix to the output file')

    h_parser = subparser.add_parser('heatmap', help='')
    h_parser.add_argument('--tsv-file', type=str, help='The TSV file that contains the extracted results')
    h_parser.add_argument('--meta-file', type=str, help='A TAB separated file that contains the sample name and conditions')
    h_parser.add_argument('--p-value', type=float, default=0.05, help='Provide a p-value cutoff (default 0.05)')
    h_parser.add_argument('--q-value', type=float, default=1.0, help='Provide a q-value cutoff (default 1.0)')
    h_parser.add_argument('--dpsi', type=float, default=0.05, help='Provide a |dPSI| cutoff (default 0.05)')
    h_parser.add_argument('--avg', type=float, default=0., help='Provide a cutoff on estimated read counts for DSA results (default 0.)')
    h_parser.add_argument('--fold-change', type=float, default=0., help='Provide a |dPSI| cutoff (default 0.)')
    h_parser.add_argument('--aggregate', action='store_true', default=False, help='Debug mode for showing more information.')
    h_parser.add_argument('--out-dir', type=str, default='./out', help='The output directory')
    h_parser.add_argument('--prefix', type=str, default='', help='Add prefix to the output file')
    h_parser.add_argument('--method', type=str, default='average',
                        help="Linkage method to use for calculating clusters. choose from 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' (default 'average')")
    h_parser.add_argument('--metric', type=str, default='braycurtis',
                        help="The distance metric to use. The distance function can be 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule' (default 'braycurtis')")

    s_parser = subparser.add_parser('sashimi', help='')
    s_parser.add_argument('--bam-list', type=str, help='A BAM files list')
    s_parser.add_argument('--tsv-file', type=str, help='The TSV file that contains the extracted results')
    s_parser.add_argument('--meta-file', type=str, help='A TAB separated file that contains the sample name and conditions')
    s_parser.add_argument('--coordinate', type=str, help='coordinate: e.g. chr1:123456-234567')
    s_parser.add_argument('--gtf', type=str, help='A GTF file')
    s_parser.add_argument('--shrink', action='store_true', default=False, help='shrink region.')
    s_parser.add_argument("--min-coverage", type=int, default=1, help='minimum coverage to ignore')
    s_parser.add_argument("--group-id", type=str, help='specify the group id')
    s_parser.add_argument("--strand", default="NONE", type=str,
        help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
    s_parser.add_argument('--out-dir', type=str, default='./out', help='The output directory')
    s_parser.add_argument('--prefix', type=str, default='', help='Add prefix to the output file')

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def run_convert_results_module(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.mntjulip_dir:
        convert_mntjulip_results(Path(args.mntjulip_dir), out_dir)
    if args.leafcutter_dir:
        convert_leafcutter_results(Path(args.leafcutter_dir), out_dir)
    if args.rmat_dir:
        convert_rmats_results(Path(args.rmat_dir), out_dir)
    if args.majiq_dir:
        convert_majiq_results(Path(args.majiq_dir), out_dir)
    if not (args.mntjulip_dir or args.leafcutter_dir or args.rmat_dir):
        raise Exception('Please specify at least a path to a program result folder!')


def run_venn_diagram_module(args):
    if not args.tsv_file_list:
        raise Exception('Please provide the list file that contains the path of the TSV result files!')

    plot_venn_diagram(Path(args.tsv_file_list), Path(args.out_dir), args.p_value, args.q_value, args.dpsi)


def run_heatmap_module(args):
    if not args.tsv_file or not args.meta_file:
        raise Exception('Please provide the list file that contains the path of the TSV result files!')
    plot_heatmap(Path(args.tsv_file), Path(args.meta_file), Path(args.out_dir), args.p_value,
                 args.q_value, args.dpsi, args.fold_change, args.avg, args.aggregate,
                 args.method, args.metric, args.prefix)


def run_sashimi_module(args):
    if args.bam_list and args.coordinate:
        sashimi_plot_with_bams(args.bam_list, args.coordinate, args.gtf, Path(args.out_dir), args.prefix, args.shrink,
                           args.strand, args.min_coverage)
    elif args.tsv_file and args.meta_file and args.gtf:
        sashimi_polt_without_bams(args.tsv_file, args.meta_file, args.gtf, args.group_id,
                                  Path(args.out_dir), args.prefix, args.shrink, args.min_coverage)
    else:
        raise Exception(textwrap.dedent('''\
    Please provide one of the two set of files,
    (1) A bam file list and a coordinate (e.g. chr1:12345-12356)
    (2) A TSV file that are generated by the 'convert-results' module, a meta file of sample names and conditions, a GTF reference genome file
    '''))


def main():
    args = get_arguments()

    if args.command == 'convert-results':
        run_convert_results_module(args)

    if args.command == 'venn-diagram':
        run_venn_diagram_module(args)

    if args.command == 'heatmap':
        run_heatmap_module(args)

    if args.command == 'sashimi':
        run_sashimi_module(args)


if __name__ == "__main__":
    main()


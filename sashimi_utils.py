import subprocess as sp
import sys, re, copy, os, codecs
from collections import OrderedDict


def sashimi_polt_without_bams(tsv_file, meta_file, gtf, group_id, out_dir, prefix, shrink, min_coverage):
    with open(tsv_file, 'r') as f:
        lines = f.readlines()

    coord, strand = None, None
    for line in lines:
        if line.startswith('#') or line.startswith('GeneName'):
            continue
        # GeneName, GroupID, FeatureElement, FeatureType, FeatureLabel, strand, p-value, q-value, dPSI,
        # ReadCount1, ReadCount2, PSI
        items = line.strip().split('\t')
        _gene_name, _group_id, label, _strand = items[0], items[1], items[4], items[5]
        if _group_id == group_id:
            strand = _strand
            chr, start, end = parse_coordinates(label)
            if coord:
                coord[1], coord[2] = min(start, coord[1]), max(end, coord[2])
            else:
                coord = [chr, start, end]
    if not coord:
        raise Exception(f"Can't find the coordinate with the provided group ID {group_id}!")
    strand = strand if strand and strand != '.' else None

    coord[1], coord[2] = coord[1] - 100, coord[2] + 100
    x, y = get_depths_from_gtf(gtf, coord, strand)

    bam_dict, overlay_dict, color_dict,  = {"+": OrderedDict()}, OrderedDict(), OrderedDict()
    label_dict, id_list = OrderedDict(), []
    if strand == '-':
        bam_dict['-'] = OrderedDict()
    if not strand:
        strand = '+'
    with open(meta_file, 'r') as f:
        for line in f:
            id, overlay_level = line.strip().split('\t')
            id_list.append(id)
            overlay_dict.setdefault(overlay_level, []).append(id)
            label_dict[overlay_level] = overlay_level
            color_dict.setdefault(overlay_level, overlay_level)
            bam_dict[strand][id] = (x, y, [], [], [], [], [])

    with open(tsv_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('#') or line.startswith('GeneName'):
            continue
        #GeneName    GroupID FeatureElement  FeatureType FeatureLabel    strand  p-value q-value dPSI    ReadCount1  ReadCount2  PSI
        items = line.strip().split('\t')
        _gene_name, _group_id, label, _strand = items[0], items[1], items[4], items[5]
        chr, _start, _end = parse_coordinates(label)
        if _group_id == group_id:
            read_counts = [int(v) for v in items[9].split(',')]
            psis = [float(v) for v in items[11].split(',')]
            chr, start, end = coord
            for i, count in enumerate(read_counts):
                # dons, accs, yd, ya, counts = [], [], [], [], []
                if count < min_coverage:
                    continue

                bam_dict[strand][id_list[i]][2].append(_start)
                bam_dict[strand][id_list[i]][3].append(_end)
                bam_dict[strand][id_list[i]][4].append( y[ _start - start - 1 ])
                bam_dict[strand][id_list[i]][5].append( y[ _end - start + 1 ])
                bam_dict[strand][id_list[i]][6].append(count)

    palette = get_preset_palette()

    # Find set of junctions to perform shrink
    intersected_introns = None
    if shrink:
        introns = (v for vs in bam_dict[strand].values() for v in zip(vs[2], vs[3]))
        intersected_introns = list(intersect_introns(introns))

    # *** PLOT *** Define plot height
    height, width, base_size = 3 * len(overlay_dict), 10, 14

    # *** PLOT *** Start R script by loading libraries, initializing variables, etc...
    R_script = setup_R_script(height, width, base_size, label_dict)
    R_script += colorize(color_dict, palette)
    R_script += make_R_lists(id_list, bam_dict[strand], overlay_dict, '', intersected_introns)

    out_format = 'png'
    out_file = out_dir / '_'.join(filter(None, [prefix, 'sashimi.png']))
    resolution = 300
    alpha = 0.5
    height = 3
    R_script += R_script_plot(out_file, out_format, resolution, '', '', height, 3, alpha)
    if os.getenv('GGSASHIMI_DEBUG') is not None:
        with open("R_script", 'w') as r:
            r.write(R_script)
    else:
        plot(R_script)


def get_depths_from_gtf(file, coord, strand):
    chr, start, end = coord
    x = [i for i in range(start, end)]
    y = [0] * (end - start)
    end = end - 1
    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            _chr, _, type, _start, _end, _, _strand, _, tags = line.strip().split("\t")
            if _chr != chr:
                continue
            _start, _end = int(_start) - 1, int(_end)
            if strand and strand != _strand:
                continue
            if type == "exon":
                if (start < _start < end or start < _end < end):
                    _start, _end = max(_start, start), min(end, _end)
                    y[_start-start: _end-start] = [1] * (_end-_start)
    return x, y

def sashimi_plot_with_bams(bams, coordinate, gtf, out_dir, prefix, shrink, strand="NONE", min_coverage=1):
    palette = get_preset_palette()

    bam_dict, overlay_dict, color_dict,  = {"+": OrderedDict()}, OrderedDict(), OrderedDict()
    label_dict, id_list = OrderedDict(), []
    if strand != "NONE":
        bam_dict["-"] = OrderedDict()

    for id, bam, overlay_level in read_bam_input(bams):
        if not os.path.isfile(bam):
            continue
        a, junctions = read_bam(bam, coordinate, strand)
        if a.keys() == ["+"] and all(map(lambda x: x == 0, list(a.values()[0]))):
            print("WARN: Sample {} has no reads in the specified area.".format(id))
            continue

        id_list.append(id)
        for _strand in a:
            bam_dict[_strand][id] = prepare_for_R(a[_strand], junctions[_strand], coordinate, min_coverage)

        overlay_dict.setdefault(overlay_level, []).append(id)
        label_dict[overlay_level] = overlay_level
        color_dict.setdefault(overlay_level, overlay_level)

    # No bam files
    if not bam_dict["+"]:
        print("ERROR: No available bam files.")
        exit(1)

    if gtf:
        transcripts, exons = read_gtf(gtf, coordinate)

    # Iterate for plus and minus strand
    for strand in bam_dict:
        # Find set of junctions to perform shrink
        intersected_introns = None
        if shrink:
            introns = (v for vs in bam_dict[strand].values() for v in zip(vs[2], vs[3]))
            intersected_introns = list(intersect_introns(introns))

        # *** PLOT *** Define plot height
        height, width, ann_height, base_size = 3 * len(overlay_dict), 10, 3, 14
        if gtf:
            height += ann_height

        # *** PLOT *** Start R script by loading libraries, initializing variables, etc...
        R_script = setup_R_script(height, width, base_size, label_dict)

        R_script += colorize(color_dict, palette)

        # *** PLOT *** Prepare annotation plot only for the first bam file
        arrow_bins = 50
        if gtf:
            # Make introns from annotation (they are shrunk if required)
            annotation = make_introns(transcripts, exons, intersected_introns)
            x = list(bam_dict[strand].values())[0][0]
            if shrink:
                x, _ = shrink_density(x, x, intersected_introns)
            R_script += gtf_for_ggplot(annotation, x[0], x[-1], arrow_bins)

        R_script += make_R_lists(id_list, bam_dict[strand], overlay_dict, '', intersected_introns)

        out_format = 'png'

        out_file = out_dir / '_'.join(filter(None, [prefix, 'sashimi.png']))
        resolution = 300
        alpha = 0.5
        height = 3
        R_script += R_script_plot(out_file, out_format, resolution, gtf, '', height, ann_height, alpha)
        if os.getenv('GGSASHIMI_DEBUG') is not None:
            with open("R_script", 'w') as r:
                r.write(R_script)
        else:
            plot(R_script)


def parse_coordinates(coord):
    coord = coord.replace(",", "")
    chr = coord.split(":")[0]
    start, end = coord.split(":")[1].split("-")
    # Convert to 0-based
    start, end = int(start) - 1, int(end)
    return chr, start, end


def count_operator(CIGAR_op, CIGAR_len, pos, start, end, array, junctions):
    # Match
    if CIGAR_op == "M":
        for i in range(pos, pos + CIGAR_len):
            if i < start or i >= end:
                continue
            ind = i - start
            array[ind] += 1

    # Insertion or Soft-clip
    if CIGAR_op == "I" or CIGAR_op == "S":
        return pos

    # Deletion
    if CIGAR_op == "D":
        pass

    # Junction
    if CIGAR_op == "N":
        don, acc = pos, pos + CIGAR_len
        if don > start and acc < end:
            junctions[(don, acc)] = junctions.setdefault((don, acc), 0) + 1

    return pos + CIGAR_len


def flip_read(s, samflag):
    if s == "NONE" or s == "SENSE":
        return 0
    if s == "ANTISENSE":
        return 1
    if s == "MATE1_SENSE":
        if int(samflag) & 64:
            return 0
        if int(samflag) & 128:
            return 1
    if s == "MATE2_SENSE":
        if int(samflag) & 64:
            return 1
        if int(samflag) & 128:
            return 0


def read_bam(file, coord, strand):
    _, start, end = parse_coordinates(coord)

    # Initialize coverage array and junction dict
    a = {"+": [0] * (end - start)}
    junctions = {"+": OrderedDict()}
    if strand != "NONE":
        a["-"] = [0] * (end - start)
        junctions["-"] = OrderedDict()

    p = sp.Popen("samtools view %s %s " % (file, coord), shell=True, stdout=sp.PIPE)
    for line in p.communicate()[0].decode('utf8').strip().split("\n"):
        if line:
            line_sp = line.strip().split("\t")
            samflag, read_start, CIGAR = line_sp[1], int(line_sp[3]), line_sp[5]

            # Ignore reads with more exotic CIGAR operators
            if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
                continue

            read_strand = ["+", "-"][flip_read(strand, samflag) ^ bool(int(samflag) & 16)]
            if strand == "NONE": read_strand = "+"

            CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
            CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]

            pos = read_start

            for n, CIGAR_op in enumerate(CIGAR_ops):
                CIGAR_len = int(CIGAR_lens[n])
                pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a[read_strand], junctions[read_strand])

    p.stdout.close()
    return a, junctions


def get_bam_path(index, path):
    if os.path.isabs(path):
        return path
    base_dir = os.path.dirname(index)
    return os.path.join(base_dir, path)


def read_bam_input(file):
    overlay = color = label = 3
    with open(file) as f:
        for line in f:
            line_sp = line.strip().split("\t")
            bam = get_bam_path(file, line_sp[1])
            overlay_level = line_sp[overlay - 1] if overlay else None
            yield line_sp[0], bam, overlay_level


def prepare_for_R(a, junctions, coord, m):
    _, start, _ = parse_coordinates(coord)

    # Convert the array index to genomic coordinates
    x = [i + start for i in range(len(a))]
    y = a

    # Arrays for R
    dons, accs, yd, ya, counts = [], [], [], [], []

    # Prepare arrays for junctions (which will be the arcs)
    for (don, acc), n in junctions.items():
        # Do not add junctions with less than defined coverage
        if n < m:
            continue

        dons.append(don)
        accs.append(acc)
        counts.append(n)

        yd.append( a[ don - start - 1 ])
        ya.append( a[ acc - start + 1 ])

    return x, y, dons, accs, yd, ya, counts


def intersect_introns(data):
    data = sorted(data)
    it = iter(data)
    a, b = next(it)
    for c, d in it:
        if b > c:  # Use `if b > c` if you want (1,2), (2,3) not to be
            # treated as intersection.
            a, b = max(a, c), min(b, d)
        else:
            yield a, b
            a, b = c, d
    yield a, b


def shrink_density(x, y, introns):
    new_x, new_y = [], []
    shift = 0
    start = 0
    # introns are already sorted by coordinates
    for a, b in introns:
        end = x.index(a) + 1
        new_x += [int(i - shift) for i in x[start:end]]
        new_y += y[start: end]
        start = x.index(b)
        l = (b - a)
        shift += l - l**0.7
    new_x += [int(i - shift) for i in x[start:]]
    new_y += y[start:]
    return new_x, new_y


def shrink_junctions(dons, accs, introns):
    new_dons, new_accs = [0] * len(dons), [0] * len(accs)
    real_introns = dict()
    shift_acc = 0
    shift_don = 0
    s = set()
    junctions = list(zip(dons, accs))
    for a, b in introns:
        l = b - a
        shift_acc += l - int(l ** 0.7)
        real_introns[a - shift_don] = a
        real_introns[b - shift_acc] = b
        for i, (don, acc) in enumerate(junctions):
            if a >= don and b <= acc:
                if (don, acc) not in s:
                    new_dons[i] = don - shift_don
                    new_accs[i] = acc - shift_acc
                else:
                    new_accs[i] = acc - shift_acc
                s.add((don,acc))
        shift_don = shift_acc
    return real_introns, new_dons, new_accs


def read_palette(f):
    palette = "#ff0000", "#00ff00", "#0000ff", "#000000"
    if f:
        with open(f) as openf:
            palette = list(line.split("\t")[0].strip() for line in openf)
    return palette


def get_preset_palette():
    # color codes from, https://sashamaps.net/docs/resources/20-colors/
    # palette = "#ff0000", "#00ff00", "#0000ff", "#000000"
    palette = '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000'
    return palette


def read_gtf(file, c):
    exons = OrderedDict()
    transcripts = OrderedDict()
    chr, start, end = parse_coordinates(c)
    end = end - 1
    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            el_chr, _, el, el_start, el_end, _, strand, _, tags = line.strip().split("\t")
            if el_chr != chr:
                continue
            d = dict(kv.strip().split(" ") for kv in tags.strip(";").split("; "))
            transcript_id = d["transcript_id"]
            el_start, el_end = int(el_start) - 1, int(el_end)
            strand = '"' + strand + '"'
            if el == "transcript":
                if (el_end > start and el_start < end):
                    transcripts[transcript_id] = max(start, el_start), min(end, el_end), strand
                continue
            if el == "exon":
                if (start < el_start < end or start < el_end < end):
                    exons.setdefault(transcript_id, []).append((max(el_start, start), min(end, el_end), strand))

    return transcripts, exons



def make_introns(transcripts, exons, intersected_introns=None):
    new_transcripts = copy.deepcopy(transcripts)
    new_exons = copy.deepcopy(exons)
    introns = OrderedDict()
    if intersected_introns:
        for tx, (tx_start, tx_end, strand) in new_transcripts.items():
            total_shift = 0
            for a, b in intersected_introns:
                l = b - a
                shift = l - int(l**0.7)
                total_shift += shift
                for i, (exon_start,exon_end,strand) in enumerate(exons.get(tx, [])):
                    new_exon_start, new_exon_end = new_exons[tx][i][:2]
                    if a < exon_start:
                        if b > exon_end:
                            if i == len(exons[tx]) - 1:
                                total_shift = total_shift - shift + (exon_start - a) * (1 - int(l**-0.3))
                            shift = (exon_start - a) * (1 - int(l**-0.3))
                            new_exon_end = new_exons[tx][i][1] - shift
                        new_exon_start = new_exons[tx][i][0] - shift
                    if b <= exon_end:
                        new_exon_end = new_exons[tx][i][1] - shift
                    new_exons[tx][i] = (new_exon_start, new_exon_end,strand)
            tx_start = min(tx_start, sorted(new_exons.get(tx, [[sys.maxsize]]))[0][0])
            new_transcripts[tx] = (tx_start, tx_end - total_shift, strand)

    for tx, (tx_start, tx_end, strand) in new_transcripts.items():
        intron_start = tx_start
        ex_end = 0
        for ex_start, ex_end, strand in sorted(new_exons.get(tx, [])):
            intron_end = ex_start
            if tx_start < ex_start:
                introns.setdefault(tx, []).append((intron_start, intron_end, strand))
            intron_start = ex_end
        if tx_end > ex_end:
            introns.setdefault(tx, []).append((intron_start, tx_end, strand))
    d = {'transcripts': new_transcripts, 'exons': new_exons, 'introns': introns}
    return d


def gtf_for_ggplot(annotation, start, end, arrow_bins):
    arrow_space = int((end - start) / arrow_bins)
    s = """

    # data table with exons
    ann_list = list(
        "exons" = data.table(),
        "introns" = data.table()
    )
    """

    if annotation["exons"]:

        s += """
        ann_list[['exons']] = data.table(
            tx = rep(c(%(tx_exons)s), c(%(n_exons)s)),
            start = c(%(exon_start)s),
            end = c(%(exon_end)s),
            strand = c(%(strand)s)
        )
        """ %({
        "tx_exons": ",".join(annotation["exons"].keys()),
        "n_exons": ",".join(map(str, map(len, annotation["exons"].values()))),
        "exon_start" : ",".join(map(str, (v[0] for vs in annotation["exons"].values() for v in vs))),
        "exon_end" : ",".join(map(str, (v[1] for vs in annotation["exons"].values() for v in vs))),
        "strand" : ",".join(map(str, (v[2] for vs in annotation["exons"].values() for v in vs))),
        })

    if annotation["introns"]:

        s += """
        ann_list[['introns']] = data.table(
            tx = rep(c(%(tx_introns)s), c(%(n_introns)s)),
            start = c(%(intron_start)s),
            end = c(%(intron_end)s),
            strand = c(%(strand)s)
        )
        # Create data table for strand arrows
        txarrows = data.table()
        introns = ann_list[['introns']]
        # Add right-pointing arrows for plus strand
        if ("+" %%in%% introns$strand && nrow(introns[strand=="+" & end-start>5, ]) > 0) {
            txarrows = rbind(
                txarrows,
                introns[strand=="+" & end-start>5, list(
                    seq(start+4,end,by=%(arrow_space)s)-1,
                    seq(start+4,end,by=%(arrow_space)s)
                    ), by=.(tx,start,end)
                ]
            )
        }
        # Add left-pointing arrows for minus strand
        if ("-" %%in%% introns$strand && nrow(introns[strand=="-" & end-start>5, ]) > 0) {
            txarrows = rbind(
                txarrows,
                introns[strand=="-" & end-start>5, list(
                    seq(start,max(start+1, end-4), by=%(arrow_space)s),
                    seq(start,max(start+1, end-4), by=%(arrow_space)s)-1
                    ), by=.(tx,start,end)
                ]
            )
        }
        """ % ({
            "tx_introns": ",".join(annotation["introns"].keys()),
            "n_introns": ",".join(map(str, map(len, annotation["introns"].values()))),
            "intron_start" : ",".join(map(str, (v[0] for vs in annotation["introns"].values() for v in vs))),
            "intron_end" : ",".join(map(str, (v[1] for vs in annotation["introns"].values() for v in vs))),
            "strand" : ",".join(map(str, (v[2] for vs in annotation["introns"].values() for v in vs))),
            "arrow_space" : arrow_space,
        })

    s += """

    gtfp = ggplot()
    if (length(ann_list[['introns']]) > 0) {
        gtfp = gtfp + geom_segment(data=ann_list[['introns']], aes(x=start, xend=end, y=tx, yend=tx), size=0.3)
        gtfp = gtfp + geom_segment(data=txarrows, aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")))
    }
    if (length(ann_list[['exons']]) > 0) {
        gtfp = gtfp + geom_segment(data=ann_list[['exons']], aes(x=start, xend=end, y=tx, yend=tx), size=5, alpha=1)
    }
    gtfp = gtfp + scale_y_discrete(expand=c(0,0.5))
    gtfp = gtfp + scale_x_continuous(expand=c(0,0.25))
    gtfp = gtfp + coord_cartesian(xlim = c(%s,%s))
    gtfp = gtfp + labs(y=NULL)
    gtfp = gtfp + theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
    """ % (start, end)

    return s


def setup_R_script(h, w, b, label_dict):
    s = """
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(data.table)
    library(gtable)

    scale_lwd = function(r) {
        lmin = 0.1
        lmax = 4
        return( r*(lmax-lmin)+lmin )
    }

    base_size = %(b)s
    height = ( %(h)s + base_size*0.352777778/67 ) * 1.02
    width = %(w)s
    theme_set(theme_bw(base_size=base_size))
    theme_update(
        plot.margin = unit(c(15,15,15,15), "pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5)
    )

    labels = list(%(labels)s)

    density_list = list()
    junction_list = list()

    """ % ({
        'h': h,
        'w': w,
        'b': b,
        'labels': ",".join(('"%s"="%s"' % (id, lab) for id, lab in label_dict.items())),
    })
    return s


def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


def mean(lst):
    return sum(lst) / len(lst)


def make_R_lists(id_list, d, overlay_dict, aggr, intersected_introns):
    s = ""
    aggr_f = {
        "mean": mean,
        "median": median,
    }
    id_list = id_list if not overlay_dict else overlay_dict.keys()
    # Iterate over ids to get bam signal and junctions
    shrinked_introns = dict()
    for k in id_list:
        shrinked_introns_k, shrinked_intronsid = dict(), dict()
        x, y, dons, accs, yd, ya, counts = [], [], [], [], [], [], []
        if not overlay_dict:
            x, y, dons, accs, yd, ya, counts = d[k]
            if intersected_introns:
                x, y = shrink_density(x, y, intersected_introns)
                shrinked_introns_k, dons, accs = shrink_junctions(dons, accs, intersected_introns)
                shrinked_introns.update(shrinked_introns_k)
        else:
            for id in overlay_dict[k]:
                xid, yid, donsid, accsid, ydid, yaid, countsid = d[id]
                if intersected_introns:
                    xid, yid = shrink_density(xid, yid, intersected_introns)
                    shrinked_intronsid, donsid, accsid = shrink_junctions(donsid, accsid, intersected_introns)
                    shrinked_introns.update(shrinked_intronsid)
                x += xid
                y += yid
                dons += donsid
                accs += accsid
                yd += ydid
                ya += yaid
                counts += countsid
            if aggr and "_j" not in aggr:
                x = d[overlay_dict[k][0]][0]
                y = list(map(aggr_f[aggr], zip(*(d[id][1] for id in overlay_dict[k]))))
                if intersected_introns:
                    x, y = shrink_density(x, y, intersected_introns)
            #dons, accs, yd, ya, counts = [], [], [], [], []
        s += """
        density_list[["%(id)s"]] = data.frame(x=c(%(x)s), y=c(%(y)s))
        junction_list[["%(id)s"]] = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))
        """ % ({
            'id': k,
            'x': ",".join(map(str, x)),
            'y': ",".join(map(str, y)),
            'dons': ",".join(map(str, dons)),
            'accs': ",".join(map(str, accs)),
            'yd': ",".join(map(str, yd)),
            'ya': ",".join(map(str, ya)),
            'counts': ",".join(map(str, counts))
        })
    if intersected_introns:
        s += """
        coord_dict = data.frame(shrinked=c(%(shrinked_introns_keys)s), real=c(%(shrinked_introns_values)s))
        intersected_introns = data.frame(real_x=c(%(intersected_introns_x)s), real_xend=c(%(intersected_introns_xend)s))
        """ % ({
            'shrinked_introns_keys': ','.join(map(str, shrinked_introns.keys())),
            'shrinked_introns_values': ','.join(map(str, shrinked_introns.values())),
            'intersected_introns_x': ','.join([str(coord[0]) for coord in intersected_introns]),
            'intersected_introns_xend': ','.join([str(coord[1]) for coord in intersected_introns])
        })
    return s


def plot(R_script):
    p = sp.Popen("R --vanilla --slave", shell=True, stdin=sp.PIPE)
    p.communicate(input=R_script.encode('utf-8'))
    p.stdin.close()
    p.wait()
    return


def colorize(color_dict, palette):
    levels = list(OrderedDict.fromkeys(color_dict.values()).keys())
    n = len(levels)
    if n > len(palette):
        palette = (palette * n)[:n]

    s = ",".join('"%s"="%s"' % (k, palette[levels.index(v)]) for k, v in color_dict.items())
    return f"color_list = list({s})\n"


def R_script_plot(file, format, resolution, gtf, aggr, height, ann_height, alpha, fix_y_scale=False):
    return """

        pdf(NULL) # just to remove the blank pdf produced by ggplotGrob

        if(packageVersion('ggplot2') >= '3.0.0'){  # fix problems with ggplot2 vs >3.0.0
            vs = 1
        } else {
            vs = 0
        }

        if(%(fix_y_scale)s) {
            maxheight = max(unlist(lapply(density_list, function(df){max(df$y)})))
            breaks_y = labeling::extended(0, maxheight, m = 4)
        }

        if(exists('coord_dict')){
            all_pos_shrinked = do.call(rbind, density_list)$x
            s2r = merge(intersected_introns, coord_dict, by.x = 'real_xend', by.y = 'real')
            s2r = merge(s2r, coord_dict, by.x = 'real_x', by.y = 'real', suffixes = c('_xend', '_x'))
            breaks_x_shrinked = labeling::extended(min(all_pos_shrinked), max(all_pos_shrinked), m = 5)
            breaks_x = c()
            out_range = c()
            for (b in breaks_x_shrinked){
                iintron = FALSE
                for (j in 1:nrow(s2r)){
                    l = s2r[j, ]
                    if(b >= l$shrinked_x && b <= l$shrinked_xend){
                        # Intersected intron
                        p = (b-l$shrinked_x)/(l$shrinked_xend - l$shrinked_x)
                        realb = round(l$real_x + p*(l$real_xend - l$real_x))
                        breaks_x = c(breaks_x, realb)
                        iintron = TRUE
                        break
                    }
                }
                if (!iintron){
                    # Exon, upstream/downstream intergenic region or intron (not intersected)
                    if(b <= min(s2r$shrinked_x)) {
                        l <- s2r[which.min(s2r$shrinked_x), ]
                        if(any(b == all_pos_shrinked)){
                            # Boundary (subtract)
                            s = l$shrinked_x - b
                            realb = l$real_x - s
                            breaks_x = c(breaks_x, realb)
                        } else {
                            out_range <- c(out_range, which(breaks_x_shrinked == b))
                        }
                    } else if (b >= max(s2r$shrinked_xend)){
                        l <- s2r[which.max(s2r$shrinked_xend), ]
                        if(any(b == all_pos_shrinked)){
                            # Boundary (sum)
                            s = b - l$shrinked_xend
                            realb = l$real_xend + s
                            breaks_x = c(breaks_x, realb)
                        } else {
                            out_range <- c(out_range, which(breaks_x_shrinked == b))
                        }
                    } else {
                        delta = b-s2r$shrinked_xend
                        delta[delta < 0] = Inf
                        l = s2r[which.min(delta), ]
                        # Internal (sum)
                        s = b - l$shrinked_xend
                        realb = l$real_xend + s
                        breaks_x = c(breaks_x, realb)
                    }
                }
            }
            if(length(out_range)) {
                breaks_x_shrinked = breaks_x_shrinked[-out_range]
            }
        }

        density_grobs = list();

        for (bam_index in 1:length(density_list)) {

            id = names(density_list)[bam_index]
            d = data.table(density_list[[id]])
            junctions = data.table(junction_list[[id]])

            # Density plot
            gp = ggplot(d) + geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=color_list[[id]], alpha=%(alpha)s)
            gp = gp + labs(y=labels[[id]])
            if(exists('coord_dict')) {
                gp = gp + scale_x_continuous(expand=c(0, 0.25), breaks = breaks_x_shrinked, labels = breaks_x)
            } else {
                gp = gp + scale_x_continuous(expand=c(0, 0.25))
            }

            if(!%(fix_y_scale)s){
                maxheight = max(d[['y']])
                breaks_y = labeling::extended(0, maxheight, m = 4)
                gp = gp + scale_y_continuous(breaks = breaks_y)
            } else {
                gp = gp + scale_y_continuous(breaks = breaks_y, limits = c(NA, maxheight))
            }

            # Aggregate junction counts
            row_i = c()
            if (nrow(junctions) >0 ) {

                junctions$jlabel = as.character(junctions$count)
                junctions = setNames(junctions[,.(max(y), max(yend),round(mean(count)),paste(jlabel,collapse=",")), keyby=.(x,xend)], names(junctions))
                if ("%(args.aggr)s" != "") {
                    junctions = setNames(junctions[,.(max(y), max(yend),round(%(args.aggr)s(count)),round(%(args.aggr)s(count))), keyby=.(x,xend)], names(junctions))
                }
                # The number of rows (unique junctions per bam) has to be calculated after aggregation
                row_i = 1:nrow(junctions)
            }

            for (i in row_i) {

                j_tot_counts = sum(junctions[['count']])

                j = as.numeric(junctions[i,1:5])

                if ("%(args.aggr)s" != "") {
                    j[3] = ifelse(length(d[x==j[1]-1,y])==0, 0, max(as.numeric(d[x==j[1]-1,y])))
                    j[4] = ifelse(length(d[x==j[2]+1,y])==0, 0, max(as.numeric(d[x==j[2]+1,y])))
                }

                # Find intron midpoint
                xmid = round(mean(j[1:2]), 1)
                ymid = max(j[3:4]) * 1.2

                # Thickness of the arch
                lwd = scale_lwd(j[5]/j_tot_counts)

                curve_par = gpar(lwd=lwd, col=color_list[[id]])

                # Arc grobs

                # Choose position of the arch (top or bottom)
                nss = i
                if (nss%%%%2 == 0) {  #bottom
                    ymid = -0.3 * maxheight
                    # Draw the arcs
                    # Left
                    curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
                    gp = gp + annotation_custom(grob = curve, j[1], xmid, 0, ymid)
                    # Right
                    curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
                    gp = gp + annotation_custom(grob = curve, xmid, j[2], 0, ymid)
                }

                if (nss%%%%2 != 0) {  #top
                    # Draw the arcs
                    # Left
                    curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
                    gp = gp + annotation_custom(grob = curve, j[1], xmid, j[3], ymid)
                    # Right
                    curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
                    gp = gp + annotation_custom(grob = curve, xmid, j[2], j[4], ymid)
                }

                # Add junction labels
                gp = gp + annotate("label", x = xmid, y = ymid, label = as.character(junctions[i,6]),
                    vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"),
                    label.size=NA, size=(base_size*0.352777778)*0.6
                )
            }

            gpGrob = ggplotGrob(gp);
            gpGrob$layout$clip[gpGrob$layout$name=="panel"] <- "off"
            if (bam_index == 1) {
                maxWidth = gpGrob$widths[2+vs] + gpGrob$widths[3+vs];    # fix problems ggplot2 vs
                maxYtextWidth = gpGrob$widths[3+vs];             # fix problems ggplot2 vs
                # Extract x axis grob (trim=F --> keep empty cells)
                xaxisGrob <- gtable_filter(gpGrob, "axis-b", trim=F)
                xaxisGrob$heights[8+vs] = gpGrob$heights[1]          # fix problems ggplot2 vs
                x.axis.height = gpGrob$heights[7+vs] + gpGrob$heights[1] # fix problems ggplot2 vs
            }

            # Remove x axis from all density plots
            kept_names = gpGrob$layout$name[gpGrob$layout$name != "axis-b"]
            gpGrob <- gtable_filter(gpGrob, paste(kept_names, sep="", collapse="|"), trim=F)

            # Find max width of y text and y label and max width of y text
            maxWidth = grid::unit.pmax(maxWidth, gpGrob$widths[2+vs] + gpGrob$widths[3+vs]); # fix problems ggplot2 vs
            maxYtextWidth = grid::unit.pmax(maxYtextWidth, gpGrob$widths[3+vs]); # fix problems ggplot2 vs
            density_grobs[[id]] = gpGrob;
        }

        # Add x axis grob after density grobs BEFORE annotation grob
        density_grobs[["xaxis"]] = xaxisGrob

        # Annotation grob
        if (%(args.gtf)s == 1) {
            gtfGrob = ggplotGrob(gtfp);
            maxWidth = grid::unit.pmax(maxWidth, gtfGrob$widths[2+vs] + gtfGrob$widths[3+vs]); # fix problems ggplot2 vs
            density_grobs[['gtf']] = gtfGrob;
        }

        # Reassign grob widths to align the plots
        for (id in names(density_grobs)) {
            density_grobs[[id]]$widths[1] <- density_grobs[[id]]$widths[1] + maxWidth - (density_grobs[[id]]$widths[2+vs] + maxYtextWidth); # fix problems ggplot2 vs
            density_grobs[[id]]$widths[3+vs] <- maxYtextWidth # fix problems ggplot2 vs
        }

        # Heights for density, x axis and annotation
        heights = unit.c(
            unit(rep(%(signal_height)s, length(density_list)), "in"),
            x.axis.height,
            unit(%(ann_height)s*%(args.gtf)s, "in")
            )

        # Arrange grobs
        argrobs = arrangeGrob(
            grobs=density_grobs,
            ncol=1,
            heights = heights,
        );

        # Save plot to file in the requested format
        if ("%(out_format)s" == "tiff"){
            # TIFF images will be lzw-compressed
            ggsave("%(out)s", plot = argrobs, device = "tiff", width = width, height = height, units = "in", dpi = %(out_resolution)s, compression = "lzw", limitsize = FALSE)
        } else {
            ggsave("%(out)s", plot = argrobs, device = "%(out_format)s", width = width, height = height, units = "in", dpi = %(out_resolution)s, limitsize = FALSE)
        }

        dev.log = dev.off()

        """ %({
            "out": file,
            "out_format": format,
            "out_resolution": resolution,
            "args.gtf": float(bool(gtf)),
            "args.aggr": aggr.rstrip("_j"),
            "signal_height": height,
            "ann_height": ann_height,
            "alpha": alpha,
            "fix_y_scale": ("TRUE" if fix_y_scale else "FALSE")
            })


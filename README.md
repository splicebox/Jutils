# Jutils

Jutils (alpha_1.0) is a visualization toolkit for alternative splicing events. Jutils directly supports vizualizing results generated by the differential splicing (DS) detection tools MntJULiP, LeafCutter, MAJIQ and rMATS, and can be easily adapted to use with other DS software.

Described in:

   Yang G, Cope L, He Z, and Florea L. (2020) Jutils: A visualization toolkit for differential alternative splicing events, *Submitted*.

```
Copyright (C) 2020-2021, and GNU GPL v3.0, by Guangyu Yang, Liliana Florea
```

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  

### <a name="table-of-contents"></a> Table of contents
- [What is Jutils?](#what-is-jutils)
- [Installation](#installation)
- [Usage](#usage)
- [Input/Output](#inputoutput)
- [Support](#support)

### <a name="what-is-jutils" /> What is Jutils?
Jutils is a visualization toolkit for alternative splicing events. It uses an intermediate tab separated file (TSV) to represent alternative splicing data such as that generated by a differential splicing prediction program. Jutils supports the vizualization of results generated by popular differential splicing (DS) prediction tools MntJULiP, LeafCutter, MAJIQ and rMATs. Additionally, users can write their own routines to convert the output of any other DS program into the unified TSV data format, which can then be visualized with Jutils.

Jutils provides three types of visualization: *heatmaps*, *sashimi plots*, and *Venn diagrams*. The user can display all records in the file, or can apply filters to select specific subsets of the data, for instance entries satisfying user defined criteria of statistical significance.

### <a name="installation"></a> Installation
MntJULiP is written in Python. To download the source code, clone the current GitHub repository:

```
git clone https://github.com/Guangyu-Yang/Jutils.git
```

#### System requirements
* Linux or Darwin  
* Python 3.7 or later

#### Prerequisites
Required Python packages: `pandas`, `numpy`, `seaborn`, `matplotlib`. The Python packages can be installed with the command:   
```
pip install --user pandas numpy seaborn matplotlib
```

### <a name="usage" />  Usage
Jutils works in two steps. *Step 1* generates the TSV file representation of the user data. *Step 2* uses the TSV file, along with other information optionally provided by the user, to generate visualizations. The basic commands are listed below, and examples as applied to the specific visualizations are given in the subsequent section.


#### Basic usage

_TSV file conversion_
```
python3 jutils.py convert-results [ --mntjulip-dir | --leafcutter-dir | --majiq-dir | rmats-dir ] <dsprogram_result_dir>
```
The command takes as input the directory containing the output from the specified DS tool, and generates a TSV file in the current(?) directory. 

_Heatmap visualization_
```
python3 jutils.py heatmap --tsv-file <tsv_file> --meta-file <meta_file> [options]
        options:
        --dpsi         cutoff for delta PSI (Percent Splice In)
        --p-value      cutoff for differential test p-value
        --q-value      cutoff for differential test q-value
        --aggregate    show results at group level (one entry per group)
        --avg          cutoff for estimated read counts of DSA results
        --fold-change  cutoff for log2(fold-change) of DSA results
        --prefix       add prefix to the output file
        --out-dir      specify the output directory
        --method       linkage method for calculating clusters
        --metric       the distance metric to use
        
```
The meta-file lists the condition for each sample, for representation in the heatmap, for instance when the data has been generated from a differential analysis. The p-value, q-value and dPSI values are those generated by the DS tool and stored in the TSV file. The format of the TSV file and that of the metadata file are shown [below](#inputoutput). Some programs including LeafCutter and MntJULiP represent introns as part of a group, and may report multiple DS introns per group. The `--aggregate` option selects one entry per group to include in the displays.

Ueser may want to visualize the differential splicing abundance (DSA) resutls generated by MntJULiP, they can apply option `--avg` to filter out introns with low estimated read counts and `--fold-change` to select the significant changed events.

`--method` and `--metric` are two options for hierarchical clustering, See [scipy.cluster.hierarchy.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage) and [scipy.spatial.distance.pdist](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist) for more information.


_Sashimi visualization_
```
python3 jutils.py sashimi --tsv-file <tsv_file> --meta-file <meta_file> [ --group-id <group_id> | --coordinate <coords> ] [options]
        where:
        --group_id     select group to visualize
        --coordinate   select genomic range for visualization
        options:
        --gtf          GTF file of gene annotations to add to the display
```
The streamlined command above will use solely information contained in the TSV and metadata files for visualization. Introns will be displayed with the values (e.g., read counts, PSI value) in multiple samples, separately by the conditions specified in the metadata file. Flanking introns will be shown at a fixed coverage level (15). Alternatively, the command below extracts the information on flanking exon coverage from the BAM files, when provided by the user, and generates traditional sashimi visualizations:

```
python3 jutils.py sashimi --bam-list <bam_list.tsv> [ --group-id <group_id> | --coordinate <coords> ] [options]
        where:
        --group_id     select group to visualize
        --coordinate   select genomic range for visualization
        options:
        --gtf          GTF file with gene annotations to add to the display
```
The format of the BAM list file is shown [below](#inputoutput).

_Venn diagram visualization_
```
python3 jutils.py venn-diagram --tsv-file-list <tsv_file_list>
```
The command above creates Venn diagram representations of gene sets from multiple TSV files, corresponding to multiple tools or comparisons. The format of the TSV file list is shown [below](#inputoutput).

#### Usage customized by program

##### LeafCutter
[LeafCutter](https://github.com/davidaknowles/leafcutter) represents alternative splicing events at the intron level, with introns having common endpoints further organized into a group. The program uses the differential splicing ratio of an intron within its group to identify DS events. The Percent Splice In (PSI) value of each intron in its group, as calculated from the raw read counts, are then used by Jutils to generate heatmaps. The user can specify criteria for inclusion, such as p-value or q-value cutoffs. Jutils produces both raw and Z-score normalized heatmaps. 

The following is an example of a script for processing results from LeafCutter:

```
# convert LeafCutter output to TSV format
result_dir='path/to/LeafCutter_results_dir'
python3 jutils.py convert-results --leafcutter-dir ${result_dir}

# heatmap: visualize DSR events
python3 jutils.py heatmap --tsv-file 'leafcutter_DSR_results.tsv' \
                          --meta-file 'meta_file.tsv' \
                          --dpsi 0.2 --p-value 0.05 --q-value 0.05                       
```

##### MntJULiP
[MntJULiP](https://github.com/splicebox/MntJULiP/) detects differentially splicing events at the level of introns, using two types of tests. The first is a differential splicing ratio (DSR) test of an intron within its group. The second is a differential test on the intron abundance level (DSA). We first describe the use of Jutils for MntJULiP(DSR), and then describe changes corresponding to the MntJULiP(DSA) test.

MntJULiP(DSR) evaluates each intron within its group ('bunch'). A 'bunch' collects all introns that share the same splice junction. PSI values for representation in the heatmap are calculated from the raw read counts in each sample. MntJULiP may report more than one DS intron per group; the `--aggregate` option can be used to display a representative value per group.

MntJULiP(DSA) evaluates each intron independently. The heatmaps represent abundance values calculated from the read counts in each sample. Instead of the dPSI value the program reports, and the TSV file stores, the change in intron abundance between the conditions compared. The user can filter records by criteria such as p-value and q-value, or the minimum fold change. The TSV file format for DSA data is described [below](#inputoutput).

Lastly, unlike all other programs, MntJULiP can perform both traditional pairwise comparisons, as well as comparisons of multiple conditions.  

Examples of scripts for all of the above scenarios follow.

_Pairwise comparison:_
```
# generate TSV file from MntJULiP output
result_dir='path/to/MntJULiP_results_dir'
python3 jutils.py convert-results --mntjulip-dir ${result_dir}

# heatmap: DSR events aggregated by groups/clusters
python3 jutils.py heatmap --tsv-file 'mntjulip_DSR_results.tsv' \
                          --meta-file 'meta_file.tsv' \
                          --dpsi 0.2 --p-value 0.05 --q-value 0.05 --aggregate

# heatmap: DSR events without aggregation
python3 jutils.py heatmap --tsv-file 'mntjulip_DSR_results.tsv' \
                          --meta-file 'meta_file.tsv' \
                          --dpsi 0.2 --p-value 0.05 --q-value 0.05

# heatmap: DSA events
python3 jutils.py heatmap --tsv-file 'mntjulip_DSA_results.tsv' \
                          --meta-file 'meta_file.tsv' \
                          --p-value 0.05 --q-value 0.05 --avg 200 --fold-change 3
```

_Multi-way comparison:_
```
generate TSV file from MntJULiP output
result_dir='path/to/MntJULiP_results_dir'
python3 jutils.py convert-results --mntjulip-dir ${result_dir}

# heatmap: DSR without aggregation
python3 jutils.py heatmap --tsv-file 'mntjulip_DSR_results.tsv' \
                          --meta-file 'meta_file_multi.tsv' \
                          --dpsi 0.2 --p-value 0.05 --q-value 0.05

# heatmap: DSA
python3 jutils.py heatmap --tsv-file 'mntjulip_DSA_results.tsv' \
                          --meta-file 'meta_file_multi.tsv' \
                          --p-value 0.05 --q-value 0.05 --avg 200
```

Sashimi visualizations use the raw read counts per sample, as before:
```
####### sashimi plot : 
python3 jutils.py convert-results --mntjulip-dir ${result_dir}

# sashimi plot with bams
python3 jutils.py sashimi --bam-list 'bam_list.tsv' \
                          --coordinate 'chr19:35744400-35745150'

# sashimi plot without bams
python3 jutils.py sashimi --meta-file 'meta_file.tsv' \
                          --tsv-file 'mntjulip_DSR_results.tsv' \
                          --group-id 'g001499' \
                          --gtf 'gencode.v27.transcripts.gtf'
```

##### rMATS
Unlike LeafCutter, MntJULiP and MAJIQ, which are intron-oriented, [rMATS](http://rnaseq-mats.sourceforge.net) reports differential splicing of canonical alternative splicing events (ASEs), including exon skipping (SE), mutually exclusive exons (MXE), alternative 5' and 3' splice sites (A5SS and A3SS) and intron retention (RI). The heatmaps represent per sample PSI values of the events calculated from the read counts. As with the other programs, intron read counts are shown on the sashimi plots.

```
# convert rMATS output to TSV format
result_dir='path/to/rMATS_results_dir'
python3 jutils.py convert-results --rmats-dir ${result_dir}

# heatmap: visualize DSR events
python3 jutils.py heatmap --tsv-file 'rmats_DSR_results.tsv' \
                          --meta-file 'meta_file.tsv' \
                          --p-value 0.05 --q-value 0.05                       
```
##### Gene set comparison
Jutils provides visualizations of gene sets generated by different comparisons or methods, and represented in TSV files, as a Venn diagram.

```
####### Venn diagram
python3 jutils.py convert-results --leafcutter-dir '/path/to/LeafCutter_results/' \
                                  --mntjulip-dir '~/path/to/MntJULiP_results/' \
                                  --rmat-dir '~/path/to/rMATS_results/'

python3 jutils.py venn-diagram --tsv-file-list tsv_file_list.txt
```
The format of the TSV file list is here.

### <a name="inputoutput" /> Input/Output

#### The unified TSV file
The TSV file contains 14 columns, describing the following attributes (for DSR data):

```
GeneName GroupID  FeatureID   FeatureType FeatureLabel   strand   p-value  q-value  dPSI  ReadCount1  ReadCount2  PSI   psi(c1)   psi(c2)
#
# LeafCutter
Mrpl15   clu_1091_NA i001  intron   chr1:4774516-4777525  .  0.788252 0.927378 -0.00266183 4,1,3,3,0,2,0,1,0,0,2,0,0,0,0,1,6,0,1,3,11,2 .  0.0444444,0.0238095,0.0285714,0.0309278,0,0.0134228,0,0.00819672,0,0,0.0194175,0,0,0,0,0.0232558,0.0451128,0,0.04,0.0191083,0.0552764,0.0135135 0.0186777433863672   0.0160159139771536
# 
# rMATS
Camk2b  .  44234   SE   chr11:5979721,5981069-5981114,5982500   -   4.84279e-13   1.71877e-09  -0.053  781,3171,1429,1317,2327,1586,2342,999,1562,1184,2082,1476,2014,2328,1568,2006,859,168,2971,2694,2703,213        75,299,213,114,489,174,270,115,183,112,180,109,157,177,104,161,65,11,243,203,199,11     0.867,0.869,0.808,0.879,0.749,0.851,0.845,0.845,0.843,0.869,0.879,0.895,0.889,0.892,0.904,0.887,0.892,0.905,0.885,0.893,0.895,0.924   .   .
```

The following slightly modified format is used for DSA type data:
```
GeneName GroupID  FeatureID  FeatureType  FeatureLabel   strand   p-value  q-value  log2FoldChange ReadCount1  ReadCount2  PSI   avg_read_counts(c1)        avg_read_counts(c2)
Xkr4   .   .  intron  chr1:3207317-3213439  -   0.383576   1   -0.22461    76,26,66,51,45,62,22,8,96,60,105,8,26,61,44,34,92,51,55,24,51,25   .    .   55.35  47.37
```

#### The meta-file
The meta-file is a TAB ('\t') separated file that lists the sample names and conditions, for example:
```
sample1  ctrl
sample2  ctrl
sample3  case
sample4  case
```
#### The bam_list file
The bam list file is a TAB ('\t') separated file containing the full path to the BAM alignment file, along with the sample name and condition, as used by the (traditional) sashimi visualization:
```
sample1   /path/to/BAM_1   ctrl
sample2   /path/to/BAM_2   case
...
```
#### The tsv_file_list
The TSV file list for Venn diagram visualizations contains the full paths of the TSV files generated for the various comparisons.
```
/path/to/TSV_for_method_1
/path/to/TSV_for_method_2
/path/to/TSV_for_method_3
```

### <a name="support" /> Support
Contact: gyang22@jhu.edu

### License information
See the file LICENSE for information on the history of this software, terms
& conditions for usage, and a DISCLAIMER OF ALL WARRANTIES.

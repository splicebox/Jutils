#! /bin/bash -l

git clone https://github.com/splicebox/Jutils.git
cd Jutils

base_dir="data"
out_dir="jutils_out"
python3 jutils.py convert-results --leafcutter-dir "${base_dir}/leafcutter" \
                                  --rmats-dir "${base_dir}/rmats" \
                                  --majiq-dir "${base_dir}/majiq" \
                                  --mntjulip-dir "${base_dir}/mntjulip" \
                                  --out-dir "${out_dir}"
result_dir="jutils_out"
# Heatmap
# DSR with aggregate
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --aggregate --prefix mntjulip_DSR

# DSR without aggregate
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR

# DSA
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSA_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --avg 50 --fold-change 2 --prefix mntjulip_DSA

# Heatmap
# DSR with aggregate, unsupervised
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --aggregate --prefix mntjulip_DSR \
                          --unsupervised

# DSR without aggregate, unsupervised
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --prefix mntjulip_DSR \
                          --unsupervised

# DSA, unsupervised
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSA_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --prefix mntjulip_DSA \
                          --unsupervised

# Venn Diagram
tsv_file="${result_dir}/tsv_file_list_1.txt"
echo "${result_dir}/mntjulip_DSR_results.tsv" > ${tsv_file}
echo "${result_dir}/leafcutter_results.tsv" >> ${tsv_file}
echo "${result_dir}/majiq_results.tsv" >> ${tsv_file}
echo "${result_dir}/rmats_ReadsOnTargetAndJunctionCounts_results.tsv" >> ${tsv_file}

python3 jutils.py venn-diagram --tsv-file-list ${tsv_file} \
                               --out-dir "${out_dir}" \
                               --p-value 0.05 --q-value 1 --dpsi 0.05 \
                               --prefix test1

# custom option for each tools
tsv_file="${result_dir}/tsv_file_list_2.txt"
echo -e "${result_dir}/mntjulip_DSR_results.tsv\tMntJULiP\t'--p-value 0.05 --q-value 1 --dpsi 0.05'" > ${tsv_file}
echo -e "${result_dir}/leafcutter_results.tsv\tLeafCutter\t'--p-value 0.05 --q-value 1 --dpsi 0.05'" >> ${tsv_file}
echo -e "${result_dir}/majiq_results.tsv\tMAJIQ\t'--p-value 0.05 --q-value 1 --dpsi 0.05'" >> ${tsv_file}
echo -e "${result_dir}/rmats_ReadsOnTargetAndJunctionCounts_results.tsv\trMATs\t'--p-value 0.05 --q-value 1 --dpsi 0.05'" >> ${tsv_file}

python3 jutils.py venn-diagram --tsv-file-list ${tsv_file} \
                               --out-dir "${out_dir}" \
                               --prefix test2



# custom option for each tools
tsv_file="${result_dir}/tsv_file_list_3.txt"
echo -e "${result_dir}/mntjulip_DSR_results.tsv\tMntJULiP1\t'--p-value 0.05 --q-value 1 --dpsi 0.05'" > ${tsv_file}
echo -e "${result_dir}/mntjulip_DSR_results.tsv\tMntJULiP2\t'--p-value 0.05 --q-value 1 --dpsi 0.1'" >> ${tsv_file}
echo -e "${result_dir}/mntjulip_DSR_results.tsv\tMntJULiP3\t'--p-value 0.05 --q-value 1 --dpsi 0.2'" >> ${tsv_file}

python3 jutils.py venn-diagram --tsv-file-list ${tsv_file} \
                               --out-dir "${out_dir}" \
                               --prefix test3

# Sashimi Plot
# MntJULiP: (used g006855 for jutils paper)
python3 jutils.py sashimi --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --group-id "g006855" \
                          --gtf "${base_dir}/gencode.vM17.annotation.clean.gtf.gz" \
                          --shrink

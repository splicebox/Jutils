
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
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --aggregate --prefix mntjulip_DSR \
                          --metric "cityblock" --method "average"

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


# Venn Diagram
/bin/ls ${result_dir}/*_results.tsv > ${result_dir}/tsv_file_list.txt
python3 jutils.py venn-diagram --tsv-file-list "${result_dir}/tsv_file_list.txt" \
                               --out-dir "${out_dir}" \
                               --p-value 0.05 --q-value 1 --dpsi 0.05

# Sashimi Plot
# MntJULiP: (used g006855 for jutils paper)
python3 jutils.py sashimi --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --group-id "g006855" \
                          --gtf "gencode.vM17.annotation.gtf" \
                          --shrink


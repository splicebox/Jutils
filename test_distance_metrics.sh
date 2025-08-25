#!/bin/bash

# Test script for distance metric support in heatmap and PCA utilities
# This script tests both existing heatmap functionality and new PCA distance metric support

set -e  # Exit on any error

# Set up directories
base_dir="$(pwd)/data"
result_dir="${base_dir}/mntjulip.jutils_testdata"
out_dir="$(pwd)/test_output"

# Create output directory
mkdir -p "${out_dir}"

echo "=== Testing Distance Metric Support in Jutils ==="
echo "Base directory: ${base_dir}"
echo "Result directory: ${result_dir}"
echo "Output directory: ${out_dir}"
echo ""

# Check if required files exist
echo "Checking for required files..."
if [[ ! -f "${result_dir}/mntjulip_DSR_results.tsv" ]]; then
    echo "Error: DSR results file not found: ${result_dir}/mntjulip_DSR_results.tsv"
    exit 1
fi

if [[ ! -f "${result_dir}/mntjulip_DSA_results.tsv" ]]; then
    echo "Error: DSA results file not found: ${result_dir}/mntjulip_DSA_results.tsv"
    exit 1
fi

if [[ ! -f "${base_dir}/mntjulip_meta_file.tsv" ]]; then
    echo "Error: Meta file not found: ${base_dir}/mntjulip_meta_file.tsv"
    exit 1
fi

echo "All required files found."
echo ""

echo "=== TESTING HEATMAP FUNCTIONALITY ==="
echo ""

echo "1. Testing DSR with aggregate (default cityblock metric)..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --aggregate --prefix mntjulip_DSR_agg
echo "✓ DSR with aggregate completed"
echo ""

echo "2. Testing DSR without aggregate (default cityblock metric)..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR
echo "✓ DSR without aggregate completed"
echo ""

echo "3. Testing DSA (default cityblock metric)..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSA_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --avg 50 --fold-change 2 --prefix mntjulip_DSA
echo "✓ DSA completed"
echo ""

echo "4. Testing DSR with aggregate, unsupervised..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --aggregate --prefix mntjulip_DSR_agg_unsup \
                          --unsupervised
echo "✓ DSR with aggregate, unsupervised completed"
echo ""

echo "5. Testing DSR without aggregate, unsupervised..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --prefix mntjulip_DSR_unsup \
                          --unsupervised
echo "✓ DSR without aggregate, unsupervised completed"
echo ""

echo "6. Testing DSA, unsupervised..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSA_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --prefix mntjulip_DSA_unsup \
                          --unsupervised
echo "✓ DSA unsupervised completed"
echo ""

echo "7. Testing heatmap with euclidean distance metric..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR_euclidean \
                          --metric euclidean
echo "✓ Heatmap with euclidean metric completed"
echo ""

echo "8. Testing heatmap with correlation distance metric..."
python3 jutils.py heatmap --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                          --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                          --out-dir "${out_dir}" \
                          --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR_correlation \
                          --metric correlation
echo "✓ Heatmap with correlation metric completed"
echo ""

echo "=== TESTING PCA FUNCTIONALITY WITH DISTANCE METRICS ==="
echo ""

echo "9. Testing PCA with default euclidean metric (standard PCA)..."
python3 jutils.py pca --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                      --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                      --out-dir "${out_dir}" \
                      --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR_pca_euclidean \
                      --metric euclidean
echo "✓ PCA with euclidean metric (standard PCA) completed"
echo ""

echo "10. Testing PCA with cityblock metric (MDS-based)..."
python3 jutils.py pca --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                      --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                      --out-dir "${out_dir}" \
                      --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR_pca_cityblock \
                      --metric cityblock
echo "✓ PCA with cityblock metric (MDS-based) completed"
echo ""

echo "11. Testing PCA with correlation metric (MDS-based)..."
python3 jutils.py pca --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                      --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                      --out-dir "${out_dir}" \
                      --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR_pca_correlation \
                      --metric correlation
echo "✓ PCA with correlation metric (MDS-based) completed"
echo ""

echo "12. Testing PCA with cosine metric (MDS-based)..."
python3 jutils.py pca --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                      --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                      --out-dir "${out_dir}" \
                      --p-value 0.05 --q-value 1 --dpsi 0.2 --prefix mntjulip_DSR_pca_cosine \
                      --metric cosine
echo "✓ PCA with cosine metric (MDS-based) completed"
echo ""

echo "13. Testing PCA with DSA data and cityblock metric..."
python3 jutils.py pca --tsv-file "${result_dir}/mntjulip_DSA_results.tsv" \
                      --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                      --out-dir "${out_dir}" \
                      --p-value 0.05 --q-value 1 --avg 50 --fold-change 2 --prefix mntjulip_DSA_pca_cityblock \
                      --metric cityblock
echo "✓ PCA with DSA data and cityblock metric completed"
echo ""

echo "14. Testing PCA unsupervised with cityblock metric..."
python3 jutils.py pca --tsv-file "${result_dir}/mntjulip_DSR_results.tsv" \
                      --meta-file "${base_dir}/mntjulip_meta_file.tsv" \
                      --out-dir "${out_dir}" \
                      --prefix mntjulip_DSR_pca_cityblock_unsup \
                      --metric cityblock \
                      --unsupervised
echo "✓ PCA unsupervised with cityblock metric completed"
echo ""

echo "=== TEST SUMMARY ==="
echo ""
echo "All tests completed successfully!"
echo ""
echo "Generated files in ${out_dir}:"
ls -la "${out_dir}"
echo ""
echo "Files with different distance metrics:"
echo "- Heatmap files with cityblock (default), euclidean, and correlation metrics"
echo "- PCA files with euclidean (standard PCA), cityblock, correlation, and cosine metrics (MDS-based)"
echo ""
echo "Key differences to observe:"
echo "1. Heatmap clustering will vary based on distance metric used"
echo "2. PCA with euclidean metric uses standard PCA (files named *pca.*)"
echo "3. PCA with other metrics uses MDS (files named *mds-{metric}.*)"
echo "4. Plot titles will show 'PCA' vs 'MDS-{metric}' to indicate the method used"
echo ""
echo "✅ Distance metric testing completed successfully!"

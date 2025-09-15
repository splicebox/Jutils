# Distance Metric Support Enhancement

This document summarizes the changes made to support city block distance (Manhattan distance) and other distance metrics in the heatmap and PCA utilities.

## Summary of Changes

### 1. Heatmap Support (Already Existing)
- The heatmap functionality **already supported** city block distance and other metrics through the `--metric` parameter
- Default metric for heatmap is `cityblock` (Manhattan distance)
- All distance metrics are passed to `sns.clustermap()` for hierarchical clustering

### 2. PCA Support (New Enhancement)
- **Added** `--metric` parameter to the PCA command
- **Enhanced** the `generate_pca()` function to support different distance metrics
- When `--metric euclidean` (default): uses standard Principal Component Analysis (PCA)
- When any other metric is specified: uses Multidimensional Scaling (MDS) with the specified distance metric

## Technical Implementation

### Files Modified:
1. `jutils.py` - Added `--metric` parameter to PCA command parser
2. `heatmap_pca_utils.py` - Enhanced `generate_pca()` function

### Key Changes in `heatmap_pca_utils.py`:
- Added imports for MDS and distance calculations:
  ```python
  from sklearn.manifold import MDS
  from scipy.spatial.distance import pdist, squareform
  ```
- Modified analysis logic:
  - For `euclidean` metric: uses standard PCA
  - For other metrics: computes distance matrix using specified metric and applies MDS
- Updated plot titles to reflect the method used (PCA vs MDS-{metric})
- Updated output filenames to include the metric type

## Usage Examples

### Heatmap with City Block Distance (Already Available)
```bash
python jutils.py heatmap --tsv-file data.tsv --meta-file meta.tsv --metric cityblock
```

### PCA with Euclidean Distance (Standard PCA)
```bash
python jutils.py pca --tsv-file data.tsv --meta-file meta.tsv --metric euclidean
```

### PCA with City Block Distance (MDS-based)
```bash
python jutils.py pca --tsv-file data.tsv --meta-file meta.tsv --metric cityblock
```

## Output File Naming
- PCA with euclidean metric: `{prefix}_pca.pc1-2.png`
- PCA with other metrics: `{prefix}_mds-{metric}.pc1-2.png`
- Heatmap files remain unchanged in naming

## Backward Compatibility
- All existing functionality is preserved
- Default behavior unchanged for heatmap (cityblock)
- Default behavior for PCA is euclidean (standard PCA)
- All existing command-line arguments work as before

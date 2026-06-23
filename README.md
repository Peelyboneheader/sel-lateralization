# SEL Lateralization

Count **Slowly Expanding Lesions (SEL)** per hemisphere (Left / Right / Middle) from a T1 MRI and an SEL mask in MNI space.

## Requirements

- **FSL** ([install](https://fsl.fmrib.ox.ac.uk/fsl/docs/install/index.html))
- **Python 3.8+**
  ```bash
  pip install nibabel nilearn scipy openpyxl matplotlib
  ```

## Usage

### Single subject

```bash
./run_sel_fsl.sh \
  --case-id SUBJ001 \
  --t1 t1.nii.gz \
  --sel sel_mask.nii.gz \
  --out ./results
```

Add `--baseline baseline_mask.nii.gz` if you have one.

### Batch (CSV)

Create a CSV with columns `CaseId,T1,SEL_mask,Lesion_baseline`:

```csv
CaseId,T1,SEL_mask,Lesion_baseline
SUBJ001,/path/to/t1_001.nii.gz,/path/to/sel_001.nii.gz,/path/to/base_001.nii.gz
SUBJ002,/path/to/t1_002.nii.gz,/path/to/sel_002.nii.gz,
```

Then run:

```bash
./run_sel_fsl.sh --manifest cases.csv --out ./results
```

## Outputs

| File | Description |
|------|-------------|
| `SEL_lateralization_results.xlsx` | Lesion counts by side per subject, with TOTAL row |
| `SEL_lesion_detail.csv` | Per-lesion MNI coordinates and side assignment |
| `<CaseId>_SEL_overlay.png` | Axial slices with lesions in red, midline in cyan |
| `SEL_lateralization_chart.png` | Bar chart of L/R/M counts |
| `registered/` | MNI-registered T1 and masks for QC |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--dof` | `12` | FLIRT degrees of freedom (6=rigid, 12=affine) |
| `--cost` | `corratio` | FLIRT cost function |
| `--ref` | MNI152_T1_1mm_brain | Registration template |

## Without FSL

```bash
python sel_lateralization.py --case-id SUBJ001 --t1 t1.nii.gz --sel sel.nii.gz --out ./results
```

Uses a built-in Python affine (less accurate than FLIRT). If data is already in MNI space, add `--no-resample`.

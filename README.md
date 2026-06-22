# SEL Lateralization

Count **Slowly Expanding Lesions (SEL)** per hemisphere — Left, Right, and Middle
(midline-crossing) — from a T1 MRI and an SEL mask, in standard MNI space.

Given a subject's T1 scan, an SEL lesion mask, and (optionally) a baseline lesion
mask, this tool registers everything to the MNI152 template, splits the brain at
the mid-sagittal plane, counts how many SEL lesions fall on each side, and writes
the results to an Excel spreadsheet. Optional overlay images and a bar chart are
produced for quick visual QC.

---

## Contents

| File | Purpose |
|------|---------|
| `run_sel_fsl.sh` | **Main entry point.** Registers with FSL FLIRT, then counts. |
| `sel_lateralization.py` | Counting + Excel/visualization engine (called by the wrapper). |
| `README.md` | This file. |

---

## Requirements

- **FSL** (provides the `flirt` command) — [install instructions](https://fsl.fmrib.ox.ac.uk/fsl/docs/install/index.html). Required for the recommended workflow.
- **Python 3.8+** with these packages:
  ```bash
  pip install nibabel nilearn scipy openpyxl matplotlib
  ```
- macOS or Linux. (On Windows, use WSL.)

---

## Quick start

```bash
# 1. Make sure FSL is active in your shell
echo $FSLDIR          # should print a path
which flirt           # should print a path

# 2. Install Python dependencies
pip install nibabel nilearn scipy openpyxl matplotlib

# 3. Make the wrapper executable (once)
chmod +x run_sel_fsl.sh

# 4. Run on one subject
./run_sel_fsl.sh \
  --case-id SUBJ001 \
  --t1 t1.nii.gz \
  --sel sel_mask.nii.gz \
  --baseline baseline_mask.nii.gz \
  --out ./results
```

Results land in `./results/` (see [Outputs](#outputs)).

> If a subject has no baseline mask, omit the `--baseline` line.

---

## Inputs

| Input | Flag | Description |
|-------|------|-------------|
| T1 MRI | `--t1` | Structural T1-weighted image (`.nii` or `.nii.gz`). |
| SEL mask | `--sel` | Binary mask of slowly expanding lesions (same space as the T1). |
| Baseline mask | `--baseline` | *(Optional)* Binary mask of baseline lesions. |
| Case ID | `--case-id` | A label for the subject (used in filenames and the spreadsheet). |

All masks must be in the **same space as their T1** (i.e. the mask overlays the
T1 correctly before registration). The tool moves the masks into MNI space using
the transform it computes from the T1.

---

## Batch processing

Create a CSV named, e.g., `cases.csv` with this exact header:

```csv
CaseId,T1,SEL_mask,Lesion_baseline
SUBJ001,/path/to/t1_001.nii.gz,/path/to/sel_001.nii.gz,/path/to/base_001.nii.gz
SUBJ002,/path/to/t1_002.nii.gz,/path/to/sel_002.nii.gz,
```

Leave the last field empty (after the comma) for subjects without a baseline mask.
Then run:

```bash
./run_sel_fsl.sh --manifest cases.csv --out ./results
```

All subjects are combined into one spreadsheet.

---

## Outputs

Written to the `--out` directory:

| File | Description |
|------|-------------|
| `SEL_lateralization_results.xlsx` | **Main result.** One row per subject: SEL counts by side, baseline count, and a TOTAL row. |
| `SEL_lesion_detail.csv` | One row per individual lesion: voxel count, MNI centroid coordinates, and assigned side. |
| `SEL_lateralization_chart.png` | Bar chart of L/R/Middle counts across subjects. |
| `<CaseId>_SEL_overlay.png` | Axial slices with SEL in red and the MNI midline drawn in cyan (visual QC). |
| `registered/<CaseId>_*_MNI.nii.gz` | The MNI-registered T1 and masks. Open in FSLeyes or 3D Slicer to verify alignment. |

### Spreadsheet columns

`CaseId`, `T1`, `SEL_mask`, `SEL_number`, `Lesion_baseline`, `SEL_left`,
`SEL_right`, `SEL_middle`. `SEL_number` is a live formula equal to
left + right + middle, so it doubles as a consistency check.

---

## How it works

1. **Register T1 → MNI.** FSL FLIRT aligns the subject's T1 to the MNI152
   template (12-parameter affine by default) and saves the transform.
2. **Move the masks.** The same transform is applied to the SEL and baseline
   masks with nearest-neighbour interpolation, so label values stay intact.
3. **Split at the midline.** In MNI space the mid-sagittal plane is the world
   x = 0 plane. Each connected lesion (26-connectivity) is assigned to a side by
   the sign of its centroid's world-x coordinate:
   - x > 0 → **Right**
   - x < 0 → **Left**
   - |x| ≤ tolerance (default 2 mm) → **Middle** (crosses midline)
4. **Report.** Counts are written to Excel; per-lesion detail to CSV; QC images
   to PNG.

The left/right convention is read from each image's affine, so it is correct
whether the data is stored in RAS (e.g. nilearn) or LAS (e.g. FSL) orientation —
no manual flipping required.

---

## Options

### Wrapper (`run_sel_fsl.sh`)

| Option | Default | Description |
|--------|---------|-------------|
| `--case-id` | — | Subject label (single-subject mode). |
| `--t1` / `--sel` / `--baseline` | — | Input images (single-subject mode). |
| `--manifest` | — | CSV for batch mode (see above). |
| `--out` | `./results` | Output directory. |
| `--dof` | `12` | FLIRT degrees of freedom (`6` = rigid, `12` = affine). |
| `--cost` | `corratio` | FLIRT cost function (e.g. `normmi` for within-modality). |
| `--ref` | `$FSLDIR/.../MNI152_T1_1mm_brain` | Registration template. |
| `--python` | `python3` | Python interpreter to use. |

### Counting engine (`sel_lateralization.py`)

Two parameters can be tuned by editing `count_lesions_by_side` in the script:
`midline_tol_mm` (default 2.0 — how close to the midline counts as "Middle") and
`min_voxels` (default 3 — drops noise specks smaller than this).

---

## Tips for good results

- **Brain-extract the T1 first** (`bet t1.nii.gz t1_brain.nii.gz`) so it matches
  the brain-extracted MNI template — this improves registration accuracy.
- **Always QC the registration.** Open `registered/<CaseId>_T1_MNI.nii.gz` with
  the SEL overlay in FSLeyes or 3D Slicer before trusting the counts. The
  `<CaseId>_SEL_overlay.png` is a fast first check.
- **Run one subject before batching** to confirm paths and orientation are right.

---

## Running without FSL (fallback)

If FSL is unavailable, `sel_lateralization.py` can self-register using a built-in
Python affine. This is less robust than FLIRT and is intended only as a fallback:

```bash
python sel_lateralization.py --case-id SUBJ001 --t1 t1.nii.gz --sel sel.nii.gz --out ./results
```

If your data is **already** in MNI space (registered by another tool), skip
registration entirely:

```bash
python sel_lateralization.py --case-id SUBJ001 --t1 t1.nii.gz --sel sel.nii.gz \
  --no-resample --out ./results
```

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---------|--------------------|
| `flirt: command not found` | FSL not active in this shell. Run `echo $FSLDIR`; if empty, source FSL (`source $FSLDIR/etc/fslconf/fsl.sh`) or reopen your terminal. |
| `Cannot work out file type` | A path points to a missing or non-NIfTI file. Check the paths in your command or CSV. |
| `Permission denied` running the wrapper | Run `chmod +x run_sel_fsl.sh` first. |
| Counts look mirrored (L/R swapped) | Confirm the mask actually overlays the T1 *before* registration; a mask in the wrong space will misregister. |
| ModuleNotFoundError | Install the Python deps: `pip install nibabel nilearn scipy openpyxl matplotlib`. |

---

## Citation / contact

If you use this in published work, please acknowledge the [maintainer]. Questions
and issues: [open a GitHub issue or contact the maintainer].

## License

[Choose a license — MIT is a common, permissive default for research code.]

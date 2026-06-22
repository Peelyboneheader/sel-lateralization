#!/usr/bin/env bash
#
# run_sel_fsl.sh — FSL FLIRT registration + SEL lateralization counting.
#
# Registers each subject's T1 to MNI152 with FLIRT, applies the SAME transform
# to the SEL mask (and optional baseline mask) with nearest-neighbour interp,
# then hands the MNI-space volumes to sel_lateralization.py (--skip-registration)
# to count Left/Right/Middle and write the Excel.
#
# Runs on macOS/Linux where FSL is installed (this is your Mac, not the sandbox).
#
# ---------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------
# Single subject:
#   ./run_sel_fsl.sh \
#       --case-id SUBJ001 \
#       --t1 t1.nii.gz \
#       --sel sel_mask.nii.gz \
#       --baseline baseline_mask.nii.gz \
#       --out ./results
#
# Batch (manifest CSV with header: CaseId,T1,SEL_mask,Lesion_baseline):
#   ./run_sel_fsl.sh --manifest cases.csv --out ./results
#
# Options:
#   --dof 12            FLIRT degrees of freedom (6=rigid, 12=affine; default 12)
#   --cost corratio     FLIRT cost function (default corratio; within-modality
#                       T1->T1 you may use normmi or normcorr)
#   --ref <path>        Override MNI reference (default: $FSLDIR/data/standard/
#                       MNI152_T1_1mm_brain.nii.gz)
#   --python <cmd>      Python interpreter (default: python3)
#   --keep-intermediate Keep the .mat and intermediate files
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DOF=12
COST="corratio"
REF=""
PYTHON="python3"
OUT="./results"
KEEP=0
MANIFEST=""
CASE_ID=""
T1=""
SEL=""
BASELINE=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/sel_lateralization.py"

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-id)   CASE_ID="$2"; shift 2;;
    --t1)        T1="$2"; shift 2;;
    --sel)       SEL="$2"; shift 2;;
    --baseline)  BASELINE="$2"; shift 2;;
    --manifest)  MANIFEST="$2"; shift 2;;
    --out)       OUT="$2"; shift 2;;
    --dof)       DOF="$2"; shift 2;;
    --cost)      COST="$2"; shift 2;;
    --ref)       REF="$2"; shift 2;;
    --python)    PYTHON="$2"; shift 2;;
    --keep-intermediate) KEEP=1; shift;;
    -h|--help)   sed -n '2,40p' "$0"; exit 0;;
    *) echo "Unknown option: $1" >&2; exit 1;;
  esac
done

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if [[ -z "${FSLDIR:-}" ]]; then
  echo "ERROR: FSLDIR is not set. Source FSL first, e.g.:" >&2
  echo "  source /usr/local/fsl/etc/fslconf/fsl.sh   (or your install path)" >&2
  exit 1
fi
command -v flirt >/dev/null 2>&1 || { echo "ERROR: 'flirt' not on PATH. Is FSL sourced?" >&2; exit 1; }

if [[ -z "$REF" ]]; then
  REF="$FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz"
fi
[[ -f "$REF" ]] || { echo "ERROR: MNI reference not found: $REF" >&2; exit 1; }
[[ -f "$PY_SCRIPT" ]] || { echo "ERROR: sel_lateralization.py not found beside this script." >&2; exit 1; }

mkdir -p "$OUT"
REG_DIR="$OUT/registered"
mkdir -p "$REG_DIR"

# ---------------------------------------------------------------------------
# Build a manifest of MNI-space files for the Python stage
# ---------------------------------------------------------------------------
MNI_MANIFEST="$OUT/_mni_manifest.csv"
echo "CaseId,T1,SEL_mask,Lesion_baseline" > "$MNI_MANIFEST"

# ---------------------------------------------------------------------------
# register_one: register T1 + masks to MNI, then APPEND the resulting row of
# MNI-space paths directly to $MNI_MANIFEST. We never capture stdout from this
# function, so FLIRT's own console output can't corrupt the manifest. All FLIRT
# stdout/stderr is redirected to this script's stderr for the same reason.
#   args: case_id  t1  sel  baseline(optional)
#   returns: 0 on success (and appends a line), nonzero on skip
# ---------------------------------------------------------------------------
register_one() {
  local cid="$1" t1="$2" sel="$3" base="${4:-}"

  [[ -f "$t1" ]]  || { echo "  [skip $cid] T1 not found: $t1" >&2; return 1; }
  [[ -f "$sel" ]] || { echo "  [skip $cid] SEL not found: $sel" >&2; return 1; }

  local t1_mni="$REG_DIR/${cid}_T1_MNI.nii.gz"
  local sel_mni="$REG_DIR/${cid}_SEL_MNI.nii.gz"
  local mat="$REG_DIR/${cid}_T1_to_MNI.mat"

  echo "  [$cid] FLIRT T1 -> MNI (dof=$DOF, cost=$COST)..." >&2
  flirt -in "$t1" -ref "$REF" -out "$t1_mni" -omat "$mat" \
        -dof "$DOF" -cost "$COST" \
        -searchrx -90 90 -searchry -90 90 -searchrz -90 90 >&2

  echo "  [$cid] Applying transform to SEL mask (nearestneighbour)..." >&2
  flirt -in "$sel" -ref "$REF" -applyxfm -init "$mat" \
        -interp nearestneighbour -out "$sel_mni" >&2

  local base_mni=""
  if [[ -n "$base" && -f "$base" ]]; then
    base_mni="$REG_DIR/${cid}_baseline_MNI.nii.gz"
    echo "  [$cid] Applying transform to baseline mask..." >&2
    flirt -in "$base" -ref "$REF" -applyxfm -init "$mat" \
          -interp nearestneighbour -out "$base_mni" >&2
  fi

  # Append the MNI-space row to the manifest (CSV).
  echo "${cid},${t1_mni},${sel_mni},${base_mni}" >> "$MNI_MANIFEST"
  return 0
}

if [[ -n "$MANIFEST" ]]; then
  [[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }
  # Skip header, read CaseId,T1,SEL_mask,Lesion_baseline. strip CR (Windows CSVs)
  # and surrounding whitespace from each field.
  tail -n +2 "$MANIFEST" | tr -d '\r' | while IFS=, read -r cid t1 sel base _rest; do
    cid="$(echo "$cid" | xargs)"; t1="$(echo "$t1" | xargs)"
    sel="$(echo "$sel" | xargs)"; base="$(echo "$base" | xargs)"
    [[ -z "$cid" ]] && continue
    register_one "$cid" "$t1" "$sel" "$base" || true
  done
else
  [[ -n "$CASE_ID" && -n "$T1" && -n "$SEL" ]] || {
    echo "ERROR: provide --manifest OR (--case-id --t1 --sel)" >&2; exit 1; }
  register_one "$CASE_ID" "$T1" "$SEL" "$BASELINE" || true
fi

# Make sure at least one case registered successfully.
if [[ "$(wc -l < "$MNI_MANIFEST")" -le 1 ]]; then
  echo "ERROR: no cases were registered. Check input paths above." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Counting + Excel (registration already done -> --skip-registration)
# ---------------------------------------------------------------------------
echo "Running lateralization counting..." >&2
"$PYTHON" "$PY_SCRIPT" --manifest "$MNI_MANIFEST" --no-resample --out "$OUT"

echo "All done. Results in: $OUT" >&2

import argparse
import math
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mzml
from tqdm import tqdm


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

PROJECT_DIR = REPO_ROOT / "sample_data" / "masscube_sample"
INPUT_CSV = PROJECT_DIR / "aligned_feature_table.csv"
SINGLE_TXT_DIR = PROJECT_DIR / "single_files"
SINGLE_CSV_DIR = PROJECT_DIR / "single_files_csv"
MZML_ROOT = PROJECT_DIR / "mzML"
BLANK_MZML_PATH = MZML_ROOT / "MB_P1-A-4_01_13240.mzML"
OUTPUT_ROOT = PROJECT_DIR / "masscube_pdet"

MZML_DIR_BY_CONC = {
    100: MZML_ROOT / "100",
    50: MZML_ROOT / "50",
    10: MZML_ROOT / "10",
}

RT_TOL = 0.3
RT_DUP_TOL = 0.3
MZ_TOL = 0.01
USE_PPM = False
MZ_TOL_PPM = 20
INT_THRESHOLD = 1000
DUP_REL_TOL = 0.02
ABS_EPS = 1e-9
WINDOW_DIFF_THRESH = 0.2
WRITE_SCANLIST_COLUMNS = True
CONVERT_SINGLE_TXT_TO_CSV = True
OVERWRITE_CONVERTED_CSV = False

MZML_PATH_CACHE = {}
SPECTRA_CACHE = {}
CSV_CACHE = {}

ATTRIBUTE_NAMES = [
    "int", "area", "width", "scancount", "snr", "smoothness", "sharpness",
    "symmetry", "density", "rank", "rtshift",
]


def resolve_path(path_value, base=REPO_ROOT):
    path = Path(path_value)
    return path if path.is_absolute() else (base / path).resolve()


def configure_paths(args):
    global PROJECT_DIR, INPUT_CSV, SINGLE_TXT_DIR, SINGLE_CSV_DIR
    global MZML_ROOT, MZML_DIR_BY_CONC, BLANK_MZML_PATH, OUTPUT_ROOT

    if args.project_dir:
        PROJECT_DIR = resolve_path(args.project_dir)

    INPUT_CSV = resolve_path(args.input_csv) if args.input_csv else PROJECT_DIR / "aligned_feature_table.csv"
    SINGLE_TXT_DIR = resolve_path(args.single_txt_dir) if args.single_txt_dir else PROJECT_DIR / "single_files"
    SINGLE_CSV_DIR = resolve_path(args.single_csv_dir) if args.single_csv_dir else PROJECT_DIR / "single_files_csv"
    MZML_ROOT = resolve_path(args.mzml_root) if args.mzml_root else PROJECT_DIR / "mzML"
    BLANK_MZML_PATH = resolve_path(args.blank_mzml) if args.blank_mzml else MZML_ROOT / "MB_P1-A-4_01_13240.mzML"
    OUTPUT_ROOT = resolve_path(args.output_root) if args.output_root else PROJECT_DIR / "masscube_pdet"

    MZML_DIR_BY_CONC = {
        100: MZML_ROOT / "100",
        50: MZML_ROOT / "50",
        10: MZML_ROOT / "10",
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert MassCube single-file TXT traces and calculate P_detection peak attributes."
    )
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--input-csv", default=None)
    parser.add_argument("--single-txt-dir", default=None)
    parser.add_argument("--single-csv-dir", default=None)
    parser.add_argument("--mzml-root", default=None)
    parser.add_argument("--blank-mzml", default=None)
    parser.add_argument("--output-root", default=None)
    parser.add_argument("--skip-txt-conversion", action="store_true")
    parser.add_argument("--overwrite-converted-csv", action="store_true")
    return parser.parse_args()


def to_float(x, default=np.nan):
    if x is None:
        return default
    if isinstance(x, str) and x.strip() == "":
        return default
    val = pd.to_numeric(pd.Series([x]), errors="coerce").iloc[0]
    if pd.isna(val):
        return default
    return float(val)


def finite(x):
    val = to_float(x)
    return not pd.isna(val) and np.isfinite(val)


def intensity_is_detected(x):
    return finite(x) and float(x) > INT_THRESHOLD


def mz_match_mask(mzs, mz_target):
    if USE_PPM:
        return np.abs((mzs - mz_target) / mz_target * 1e6) <= MZ_TOL_PPM
    return np.abs(mzs - mz_target) <= MZ_TOL


def is_intensity_match(scan_int, aligned_int):
    if not finite(scan_int) or not finite(aligned_int) or float(aligned_int) <= 0:
        return False
    rel_diff = abs(float(scan_int) - float(aligned_int)) / float(aligned_int)
    return rel_diff <= DUP_REL_TOL or abs(float(scan_int) - float(aligned_int)) <= ABS_EPS


def assign_p_detection(detection_count, num_reps):
    return 0 if num_reps <= 0 else int(round(detection_count * 100.0 / float(num_reps)))


def safe_trapz(y, x):
    if len(y) < 2:
        return 0.0
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def safe_sort_for_feature_id(df):
    if "feature_ID" not in df.columns:
        return df
    temp = df.copy()
    temp["_feature_sort"] = pd.to_numeric(temp["feature_ID"], errors="coerce")
    if temp["_feature_sort"].notna().any():
        temp = temp.sort_values(["_feature_sort", "feature_ID"], kind="mergesort")
    else:
        temp = temp.sort_values("feature_ID", kind="mergesort")
    return temp.drop(columns=["_feature_sort"], errors="ignore")


def strip_internal_columns(df):
    return df.drop(columns=[c for c in df.columns if str(c).startswith("_")], errors="ignore")


def convert_single_txt_files_to_csv(txt_dir, csv_dir, overwrite=False):
    txt_dir = Path(txt_dir)
    csv_dir = Path(csv_dir)
    if not txt_dir.exists():
        print(f"Single-file TXT directory not found: {txt_dir}")
        return {"converted": 0, "skipped": 0}

    csv_dir.mkdir(parents=True, exist_ok=True)
    txt_files = sorted(p for p in txt_dir.iterdir() if p.is_file() and p.suffix.lower() == ".txt")
    converted = 0
    skipped = 0

    print(f"Converting TXT traces from {txt_dir} to {csv_dir}")
    for txt_path in tqdm(txt_files, desc="txt -> csv", ncols=90):
        out_path = csv_dir / f"{txt_path.stem}.csv"
        if out_path.exists() and not overwrite:
            skipped += 1
            continue
        df = pd.read_csv(txt_path, sep="\t", encoding="latin1", engine="python")
        if df.shape[1] == 1:
            df = pd.read_csv(txt_path, sep=None, encoding="latin1", engine="python")
        df.to_csv(out_path, index=False, encoding="latin1")
        converted += 1

    print(f"TXT conversion summary: {converted} converted, {skipped} skipped.")
    return {"converted": converted, "skipped": skipped}


def find_first_existing_col(df, candidates):
    lower_map = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        key = cand.lower()
        if key in lower_map:
            return lower_map[key]
    return None


def discover_masscube_samples(df):
    pattern = re.compile(r"^(?P<conc>\d+)-(?P<rep>\d+)_")
    groups = defaultdict(list)

    for col in df.columns:
        col_str = str(col)
        if col_str.startswith("MB-"):
            continue
        match = pattern.match(col_str)
        if not match:
            continue
        conc = int(match.group("conc"))
        rep = int(match.group("rep"))
        groups[conc].append({
            "conc": conc,
            "rep": rep,
            "sample_name": col_str,
            "sample_base": col_str,
            "intensity_col": col_str,
        })

    for conc in groups:
        groups[conc] = sorted(groups[conc], key=lambda s: s["rep"])
    return dict(groups)


def prepare_feature_table(input_csv):
    df = pd.read_csv(input_csv, dtype=str, encoding="latin1")
    rt_col = find_first_existing_col(df, ["RT", "rt", "retention_time"])
    mz_col = find_first_existing_col(df, ["m/z", "mz", "MZ"])
    fid_col = find_first_existing_col(df, ["feature_ID", "feature_id", "id"])
    group_col = find_first_existing_col(df, ["group_ID", "group_id"])

    if rt_col is None or mz_col is None:
        raise ValueError("Aligned table must contain RT and m/z columns.")

    df["_aligned_rt"] = pd.to_numeric(df[rt_col], errors="coerce")
    df["_aligned_mz"] = pd.to_numeric(df[mz_col], errors="coerce")

    if fid_col is None:
        df["_feature_id"] = np.arange(1, len(df) + 1).astype(str)
        fid_col = "_feature_id"

    return df, {
        "rt_col": rt_col,
        "mz_col": mz_col,
        "fid_col": fid_col,
        "group_col": group_col,
    }


def normalize_stem(name):
    stem = os.path.splitext(os.path.basename(name))[0]
    if stem.lower().endswith(".d"):
        stem = stem[:-2]
    return stem.lower()


def find_mzml_path_for_sample(sample, mzml_dir):
    key = (sample["sample_base"].lower(), str(mzml_dir))
    if key in MZML_PATH_CACHE:
        return MZML_PATH_CACHE[key]

    mzml_dir = Path(mzml_dir)
    if not mzml_dir.is_dir():
        raise FileNotFoundError(f"mzML directory does not exist: {mzml_dir}")

    files = sorted(p for p in mzml_dir.iterdir() if p.is_file() and p.suffix.lower() == ".mzml")
    target_norm = sample["sample_base"].lower()

    exact_matches = [p for p in files if normalize_stem(p.name) == target_norm]
    if len(exact_matches) == 1:
        MZML_PATH_CACHE[key] = exact_matches[0]
        return exact_matches[0]

    prefix = f"{sample['conc']}-{sample['rep']}_".lower()
    prefix_matches = [p for p in files if p.name.lower().startswith(prefix)]
    if len(prefix_matches) == 1:
        MZML_PATH_CACHE[key] = prefix_matches[0]
        return prefix_matches[0]

    raise FileNotFoundError(
        f"Could not uniquely match mzML for {sample['sample_name']} in {mzml_dir}. "
        f"Exact matches: {len(exact_matches)}; prefix matches: {len(prefix_matches)}."
    )


def spectrum_rt_minutes(spec):
    scan_list = spec.get("scanList", {}).get("scan", [])
    if not scan_list:
        return np.nan
    scan = scan_list[0]
    rt = to_float(scan.get("scan start time"))
    if pd.isna(rt):
        return np.nan
    unit_name = str(scan.get("scan start time unitName", scan.get("unitName", "minute"))).lower()
    if "second" in unit_name:
        return rt / 60.0
    return rt


def load_mzml_file(path):
    path = Path(path)
    if path in SPECTRA_CACHE:
        return SPECTRA_CACHE[path]

    spectra = []
    with mzml.read(str(path)) as reader:
        for scan_idx, spec in enumerate(reader):
            if spec.get("ms level") != 1:
                continue
            rt = spectrum_rt_minutes(spec)
            if pd.isna(rt):
                continue
            if "m/z array" not in spec or "intensity array" not in spec:
                continue
            mzs = np.asarray(spec["m/z array"], dtype=np.float32)
            ints = np.asarray(spec["intensity array"], dtype=np.float32)
            if mzs.size == 0 or ints.size == 0:
                continue
            spectra.append({"rt": float(rt), "mzs": mzs, "ints": ints, "scan_idx": int(scan_idx)})

    spectra = sorted(spectra, key=lambda s: s["rt"])
    SPECTRA_CACHE[path] = spectra
    return spectra


def load_mzml_spectra(sample, mzml_dir):
    path = find_mzml_path_for_sample(sample, mzml_dir)
    return load_mzml_file(path)


def parse_peak_shape_window(peak_shape_value):
    if peak_shape_value is None or pd.isna(peak_shape_value):
        return np.nan, np.nan
    tokens = [t for t in str(peak_shape_value).split("|") if ";" in t]
    if not tokens:
        return np.nan, np.nan
    first = tokens[0].split(";")
    last = tokens[-1].split(";")
    if len(first) < 2 or len(last) < 2:
        return np.nan, np.nan
    return to_float(first[0]), to_float(last[0])


def detect_area_column(df):
    candidates = ["peak_area", "area", "Area", "peak area", "Peak area", "peak_area_raw", "peak_area_smoothed"]
    for col in candidates:
        if col in df.columns:
            return col
    lower_map = {str(c).lower(): c for c in df.columns}
    for key in ["peak_area", "area", "peak area"]:
        if key in lower_map:
            return lower_map[key]
    return None


def load_sample_csv(sample_base_name):
    if sample_base_name in CSV_CACHE:
        return CSV_CACHE[sample_base_name]

    fpath = SINGLE_CSV_DIR / f"{sample_base_name}.csv"
    if not fpath.exists():
        CSV_CACHE[sample_base_name] = None
        return None

    df = pd.read_csv(fpath, encoding="latin1")
    for col in ["RT", "m/z", "RT_start", "RT_end", "peak_height", "peak_shape"]:
        if col not in df.columns:
            df[col] = "" if col == "peak_shape" else np.nan

    for col in ["RT", "m/z", "RT_start", "RT_end", "peak_height"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["peak_shape"] = df["peak_shape"].astype(str)

    area_col = detect_area_column(df)
    df["_area"] = pd.to_numeric(df[area_col], errors="coerce") if area_col else np.nan

    windows = df["peak_shape"].map(parse_peak_shape_window)
    df["PS_RT_start"] = [w[0] for w in windows]
    df["PS_RT_end"] = [w[1] for w in windows]

    CSV_CACHE[sample_base_name] = df
    return df


def find_rows_in_sample_csv(sample_df, rt_val, mz_val):
    if sample_df is None or sample_df.empty:
        return pd.DataFrame()
    rt_mask = sample_df["RT"].notna() & sample_df["RT"].between(rt_val - RT_TOL, rt_val + RT_TOL)
    if USE_PPM:
        ppm_diffs = np.abs((sample_df["m/z"] - mz_val) / mz_val * 1e6)
        mz_mask = sample_df["m/z"].notna() & (ppm_diffs <= MZ_TOL_PPM)
    else:
        mz_mask = sample_df["m/z"].notna() & sample_df["m/z"].between(mz_val - MZ_TOL, mz_val + MZ_TOL)
    return sample_df.loc[rt_mask & mz_mask].copy()


def pick_best_candidate(candidate_rows, aligned_rt, aligned_intensity=None):
    if candidate_rows is None or candidate_rows.empty:
        return None

    if aligned_intensity is not None and intensity_is_detected(aligned_intensity):
        candidate_rows = candidate_rows.copy()
        candidate_rows["_int_match"] = candidate_rows["peak_height"].map(
            lambda x: is_intensity_match(x, aligned_intensity)
        )
        matches = candidate_rows.loc[candidate_rows["_int_match"]].copy()
        if not matches.empty:
            matches["_rt_diff"] = np.abs(matches["RT"] - aligned_rt)
            return matches.sort_values("_rt_diff").iloc[0]

    candidate_rows = candidate_rows.copy()
    candidate_rows["_rt_diff"] = np.abs(candidate_rows["RT"] - aligned_rt)
    return candidate_rows.sort_values("_rt_diff").iloc[0]


def get_window_from_sample_csv(row, sample):
    aligned_rt = float(row["_aligned_rt"])
    aligned_mz = float(row["_aligned_mz"])
    aligned_int = to_float(row.get(sample["intensity_col"]))
    sample_df = load_sample_csv(sample["sample_base"])
    candidates = find_rows_in_sample_csv(sample_df, aligned_rt, aligned_mz)
    chosen = pick_best_candidate(candidates, aligned_rt, aligned_intensity=aligned_int)

    if chosen is None:
        return aligned_rt - RT_TOL, aligned_rt + RT_TOL, "aligned_rt_tol", np.nan

    ps_start = to_float(chosen.get("PS_RT_start"))
    ps_end = to_float(chosen.get("PS_RT_end"))
    if finite(ps_start) and finite(ps_end) and ps_end > ps_start:
        return float(ps_start), float(ps_end), "peak_shape", to_float(chosen.get("_area"))

    rt_start = to_float(chosen.get("RT_start"))
    rt_end = to_float(chosen.get("RT_end"))
    if finite(rt_start) and finite(rt_end) and rt_end > rt_start:
        return float(rt_start), float(rt_end), "RT_start_end", to_float(chosen.get("_area"))

    chosen_rt = to_float(chosen.get("RT"), default=aligned_rt)
    return float(chosen_rt - RT_TOL), float(chosen_rt + RT_TOL), "sample_rt_tol", to_float(chosen.get("_area"))


def expand_window_if_peakshape_edge(row, sample, rt_start, rt_end):
    aligned_rt = float(row["_aligned_rt"])
    aligned_mz = float(row["_aligned_mz"])
    aligned_int = to_float(row.get(sample["intensity_col"]))
    sample_df = load_sample_csv(sample["sample_base"])
    candidates = find_rows_in_sample_csv(sample_df, aligned_rt, aligned_mz)
    chosen = pick_best_candidate(candidates, aligned_rt, aligned_intensity=aligned_int)

    if chosen is None:
        return rt_start, rt_end

    tokens = [t for t in str(chosen.get("peak_shape", "")).split("|") if ";" in t]
    if not tokens:
        return rt_start, rt_end

    first = tokens[0].split(";")
    last = tokens[-1].split(";")
    first_int = to_float(first[1] if len(first) > 1 else np.nan)
    last_int = to_float(last[1] if len(last) > 1 else np.nan)
    peak_height = to_float(chosen.get("peak_height"), default=1.0)

    if pd.isna(first_int) or pd.isna(last_int) or peak_height <= 0:
        return rt_start, rt_end
    if last_int > peak_height * 0.1 and aligned_rt > rt_end - 0.1:
        rt_end = min(rt_end + RT_TOL, aligned_rt + RT_TOL)
    if first_int > peak_height * 0.1 and aligned_rt < rt_start + 0.1:
        rt_start = max(rt_start - RT_TOL, aligned_rt - RT_TOL)
    return rt_start, rt_end


def windows_pairwise_inconsistent(windows, thresh=WINDOW_DIFF_THRESH):
    valid = [(s, e) for s, e in windows if finite(s) and finite(e)]
    for i in range(len(valid)):
        for j in range(i + 1, len(valid)):
            if abs(valid[i][0] - valid[j][0]) > thresh or abs(valid[i][1] - valid[j][1]) > thresh:
                return True
    return False


def compute_intersection_window(windows):
    valid = [(float(s), float(e)) for s, e in windows if finite(s) and finite(e)]
    if not valid:
        return None, None
    start = max(s for s, _ in valid)
    end = min(e for _, e in valid)
    return (start, end) if end > start else (None, None)


def harmonize_windows(per_sample_windows):
    out = dict(per_sample_windows)
    good = [(rep, float(s), float(e)) for rep, (s, e, src) in out.items() if finite(s) and finite(e) and src in {"peak_shape", "RT_start_end"}]

    if len(good) >= 2:
        mean_start = float(np.mean([x[1] for x in good]))
        mean_end = float(np.mean([x[2] for x in good]))
        for rep, (s, e, src) in list(out.items()):
            if finite(s) and finite(e) and (float(s) < mean_start - WINDOW_DIFF_THRESH or float(e) > mean_end + WINDOW_DIFF_THRESH):
                out[rep] = (mean_start, mean_end, "trimmed_to_mean")
        return out

    windows = [(s, e) for s, e, _ in out.values()]
    if not windows_pairwise_inconsistent(windows):
        return out

    start, end = compute_intersection_window(windows)
    if start is None or end is None:
        valid = [(float(s), float(e)) for s, e in windows if finite(s) and finite(e)]
        if valid:
            start = min(s for s, _ in valid)
            end = max(e for _, e in valid)

    if start is not None and end is not None and end > start:
        for rep in out:
            out[rep] = (start, end, "harmonized")
    return out


def clip_window_to_mzml(rt_start, rt_end, aligned_rt, spectra_list, source):
    if not spectra_list:
        return rt_start, rt_end, source
    mzml_rt_min = spectra_list[0]["rt"]
    mzml_rt_max = spectra_list[-1]["rt"]
    if rt_end < mzml_rt_min or rt_start > mzml_rt_max:
        rt_start = max(mzml_rt_min, aligned_rt - RT_TOL)
        rt_end = min(mzml_rt_max, aligned_rt + RT_TOL)
        source = f"{source}+mzml_clip"
    return rt_start, rt_end, source


def extract_scans_within_rt_window(spectra_list, rt_start, rt_end, mz_target, exclude_scan_idxs=None):
    if not spectra_list or rt_start is None or rt_end is None or rt_end < rt_start:
        return [], [], []

    rts = np.array([s["rt"] for s in spectra_list], dtype=np.float64)
    start_idx = np.searchsorted(rts, rt_start, side="left")
    end_idx = np.searchsorted(rts, rt_end, side="right")
    exclude_scan_idxs = exclude_scan_idxs or set()
    all_scans = []

    for spec in spectra_list[start_idx:end_idx]:
        scan_idx = spec["scan_idx"]
        if scan_idx in exclude_scan_idxs:
            continue
        mask = mz_match_mask(spec["mzs"], mz_target)
        if not np.any(mask):
            continue
        all_scans.append({
            "scan_idx": scan_idx,
            "rt": float(spec["rt"]),
            "max_int": float(np.max(spec["ints"][mask])),
        })

    lowint_scans = [s for s in all_scans if not intensity_is_detected(s["max_int"])]
    highint_scans = [s for s in all_scans if intensity_is_detected(s["max_int"])]
    return all_scans, lowint_scans, highint_scans


def find_duplicate_scan_idxs(aligned_df, current_idx, sample, aligned_rt, aligned_mz, scans_to_check):
    rep_col = sample["intensity_col"]
    if rep_col not in aligned_df.columns:
        return set(), []

    rt_mask = aligned_df["_aligned_rt"].between(aligned_rt - RT_DUP_TOL, aligned_rt + RT_DUP_TOL)
    if USE_PPM:
        ppm_diffs = np.abs((aligned_df["_aligned_mz"] - aligned_mz) / aligned_mz * 1e6)
        mz_mask = ppm_diffs <= MZ_TOL_PPM
    else:
        mz_mask = aligned_df["_aligned_mz"].between(aligned_mz - MZ_TOL, aligned_mz + MZ_TOL)

    candidates = aligned_df.loc[rt_mask & mz_mask].copy()
    if current_idx in candidates.index:
        candidates = candidates.drop(index=current_idx)
    if candidates.empty:
        return set(), []

    candidates["_candidate_int"] = pd.to_numeric(candidates[rep_col], errors="coerce").fillna(0)
    candidates = candidates.loc[candidates["_candidate_int"] > 0]
    duplicate_scan_idxs = set()
    duplicate_notes = []

    for scan in scans_to_check:
        for cand_idx, cand in candidates.iterrows():
            cand_int = float(cand["_candidate_int"])
            if is_intensity_match(scan["max_int"], cand_int):
                duplicate_scan_idxs.add(scan["scan_idx"])
                duplicate_notes.append({
                    "scan_idx": scan["scan_idx"],
                    "scan_rt": scan["rt"],
                    "scan_int": scan["max_int"],
                    "duplicate_feature_index": cand_idx,
                    "duplicate_int": cand_int,
                })
                break
    return duplicate_scan_idxs, duplicate_notes


def scanlist_string(scans):
    if not scans:
        return ""
    scans = sorted(scans, key=lambda x: x["rt"])
    return "; ".join(f"{s['scan_idx'] + 1}|{s['rt']:.4f}|{int(round(s['max_int']))}" for s in scans)


def scanid_string(scans):
    if not scans:
        return ""
    scans = sorted(scans, key=lambda x: x["rt"])
    return ", ".join(str(s["scan_idx"] + 1) for s in scans)


def interp_rt_at_height(rts, ints, target_height):
    order = np.argsort(ints)
    ints_sorted = ints[order]
    rts_sorted = rts[order]
    unique_ints, unique_idx = np.unique(ints_sorted, return_index=True)
    unique_rts = rts_sorted[unique_idx]
    if len(unique_ints) < 2:
        return np.nan
    return float(np.interp(target_height, unique_ints, unique_rts))


def calculate_shape_metrics_from_scans(scans):
    if not scans:
        return 0.0, 0.0, 1.0

    scans = sorted(scans, key=lambda x: x["rt"])
    all_rts = np.array([float(s["rt"]) for s in scans], dtype=float)
    all_ints = np.array([float(s["max_int"]) for s in scans], dtype=float)
    n_total = len(all_ints)

    if n_total < 3:
        return 0.0, 0.0, 1.0

    apex_idx = int(np.argmax(all_ints))
    apex_val = float(all_ints[apex_idx])
    apex_rt = float(all_rts[apex_idx])
    if apex_val <= 0:
        return 0.0, 0.0, 1.0

    target_height = 0.10 * apex_val
    rt_left = interp_rt_at_height(all_rts[:apex_idx + 1], all_ints[:apex_idx + 1], target_height)
    rt_right = interp_rt_at_height(all_rts[apex_idx:], all_ints[apex_idx:], target_height)
    symmetry_score = 1.0
    if finite(rt_left) and finite(rt_right):
        a = apex_rt - rt_left
        b = rt_right - apex_rt
        if a > 0:
            symmetry_score = b / a

    sigma = max(1.0, n_total / 6.0)
    deltas = np.diff(all_ints)
    score_numerator = 0.0
    score_denominator = 0.0

    for i, d in enumerate(deltas):
        dist = abs(i - apex_idx)
        w_dist = np.exp(-(dist ** 2) / (2 * sigma ** 2))
        w_int = np.clip(max(all_ints[i], all_ints[i + 1]) / apex_val, 0.1, 1.0)
        abs_d = abs(d) * w_dist * w_int
        score_denominator += abs_d
        if (i < apex_idx and d > 0) or (i >= apex_idx and d < 0):
            score_numerator += abs_d

    smoothness_score = score_numerator / score_denominator if score_denominator > 0 else 0.0
    sqrt_apex = np.sqrt(apex_val)
    sharp_vals = [abs(apex_val - all_ints[i]) / (abs(apex_idx - i) * sqrt_apex) for i in range(n_total) if i != apex_idx]
    sharpness_score = float(np.max(sharp_vals)) if sharp_vals else 0.0
    return float(smoothness_score), float(sharpness_score), float(symmetry_score)


def calculate_trace_attributes(scans, rt_start, rt_end):
    width = float(rt_end - rt_start) if finite(rt_start) and finite(rt_end) and rt_end > rt_start else np.nan
    if not scans:
        return {
            "width": width,
            "scancount": 0,
            "apex_rt": 0.0,
            "raw_max": 0.0,
            "raw_area": 0.0,
            "smoothness": 0.0,
            "sharpness": 0.0,
            "symmetry": 1.0,
        }

    scans = sorted(scans, key=lambda x: x["rt"])
    rts = np.array([s["rt"] for s in scans], dtype=float)
    ints = np.array([s["max_int"] for s in scans], dtype=float)
    apex_pos = int(np.argmax(ints))
    smoothness, sharpness, symmetry = calculate_shape_metrics_from_scans(scans)

    return {
        "width": width,
        "scancount": int(len(scans)),
        "apex_rt": float(rts[apex_pos]),
        "raw_max": float(ints[apex_pos]),
        "raw_area": safe_trapz(ints, rts),
        "smoothness": smoothness,
        "sharpness": sharpness,
        "symmetry": symmetry,
    }


def load_blank_spectra_for_noise(blank_path):
    blank_path = Path(blank_path)
    if not blank_path.exists():
        print(f"Blank mzML not found: {blank_path}. S/N will use noise = 1.0.")
        return []

    blank_data = []
    with mzml.read(str(blank_path)) as reader:
        for spec in reader:
            if spec.get("ms level") != 1:
                continue
            if "m/z array" not in spec or "intensity array" not in spec:
                continue
            mzs = np.asarray(spec["m/z array"], dtype=float)
            ints = np.asarray(spec["intensity array"], dtype=float)
            if mzs.size == 0 or ints.size == 0:
                continue
            blank_data.append({"mzs": mzs, "ints": ints})

    print(f"Loaded {len(blank_data)} blank MS1 scans for S/N.")
    return blank_data


def calculate_noise_map(feature_records, blank_data):
    noise_map = {}
    if not blank_data:
        for rec in feature_records:
            noise_map[rec["feature_ID"]] = 1.0
        return noise_map

    n_nonzero = 0
    n_defaulted = 0

    for rec in tqdm(feature_records, desc="blank noise", ncols=90):
        fid = rec["feature_ID"]
        target_mz = float(rec["m/z"])
        points = []
        for spec in blank_data:
            mask = mz_match_mask(spec["mzs"], target_mz)
            points.append(float(np.max(spec["ints"][mask])) if np.any(mask) else 0.0)
        noise = float(np.std(points)) if points else 0.0
        if noise > 0:
            noise_map[fid] = noise
            n_nonzero += 1
        else:
            noise_map[fid] = 1.0
            n_defaulted += 1

    print(f"Blank noise summary: {n_nonzero} features had nonzero SD; {n_defaulted} defaulted to 1.0.")
    return noise_map


def add_density_rank_rtshift(records, reps):
    if not records:
        return records

    aligned_rts = np.array([float(r["RT"]) for r in records], dtype=float)
    for rep in reps:
        apex_rts = np.array([to_float(r.get(f"apex_rt_{rep}"), default=0.0) for r in records], dtype=float)
        context_ints = np.array([to_float(r.get(f"context_int_{rep}"), default=0.0) for r in records], dtype=float)

        for i, rec in enumerate(records):
            curr_rt = apex_rts[i]
            curr_int = context_ints[i]
            if curr_int <= 0 or curr_rt <= 0:
                rec[f"density_{rep}"] = 0
                rec[f"rank_{rep}"] = 0
                rec[f"rtshift_{rep}"] = 0.0
                continue
            mask = (apex_rts >= curr_rt - RT_TOL) & (apex_rts <= curr_rt + RT_TOL) & (context_ints > 0)
            local_ints = context_ints[mask]
            rec[f"density_{rep}"] = int(len(local_ints))
            rec[f"rank_{rep}"] = int((local_ints > curr_int).sum() + 1)
            rec[f"rtshift_{rep}"] = float(curr_rt - aligned_rts[i])
    return records


def build_clean_attribute_table(full_df, reps):
    meta_cols = [c for c in ["group_ID", "feature_ID", "RT", "m/z", "p_detection"] if c in full_df.columns]
    ordered_cols = meta_cols.copy()
    for attr in ATTRIBUTE_NAMES:
        ordered_cols.extend([f"{attr}_{r}" for r in reps if f"{attr}_{r}" in full_df.columns])
    return full_df[[c for c in ordered_cols if c in full_df.columns]].copy()


def build_extraction_analysis_table(full_df, reps):
    meta_cols = [c for c in [
        "group_ID", "feature_ID", "RT", "m/z", "p_detection",
        "detection_count", "num_replicates", "not_detected_rep_original", "blank_noise_sd",
    ] if c in full_df.columns]
    prefixes = [
        "original_int", "original_detected", "int", "area", "scancount", "apex_rt",
        "context_int", "raw_max", "raw_area", "rt_start", "rt_end", "window_source",
        "int_source", "area_source", "excluded_low", "excluded_dup", "scanid", "scanlist",
    ]
    cols = meta_cols.copy()
    for prefix in prefixes:
        cols.extend([f"{prefix}_{r}" for r in reps if f"{prefix}_{r}" in full_df.columns])
    return full_df[[c for c in cols if c in full_df.columns]].copy()


def write_attribute_files(clean_df, reps, output_dir):
    attr_dir = Path(output_dir) / "attributes"
    attr_dir.mkdir(parents=True, exist_ok=True)
    meta_cols = [c for c in ["group_ID", "feature_ID", "RT", "m/z", "p_detection"] if c in clean_df.columns]
    for attr in ATTRIBUTE_NAMES:
        cols = meta_cols + [f"{attr}_{r}" for r in reps if f"{attr}_{r}" in clean_df.columns]
        if len(cols) > len(meta_cols):
            clean_df[cols].to_csv(attr_dir / f"{attr}.csv", index=False)


def write_summary(clean_df, output_dir):
    summary_path = Path(output_dir) / "masscube_pdet_summary.csv"
    if clean_df.empty:
        pd.DataFrame(columns=["p_detection", "feature_count"]).to_csv(summary_path, index=False)
        return
    summary = clean_df.groupby("p_detection", dropna=False).size().reset_index(name="feature_count").sort_values("p_detection")
    summary.to_csv(summary_path, index=False)


def process_concentration(df, cols, conc, samples, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    strip_internal_columns(df).to_csv(output_dir / "original_table.csv", index=False)

    mzml_dir = MZML_DIR_BY_CONC.get(conc)
    if mzml_dir is None:
        raise ValueError(f"No mzML directory configured for concentration {conc}.")

    reps = [s["rep"] for s in samples]
    n_reps = len(samples)
    fid_col = cols["fid_col"]
    group_col = cols.get("group_col")

    print(f"\nProcessing concentration {conc} with {n_reps} samples")
    print(f"Aligned table: {INPUT_CSV}")
    print(f"Sample CSV directory: {SINGLE_CSV_DIR}")
    print(f"mzML directory: {mzml_dir}")
    print(f"Output directory: {output_dir}")

    spectra_by_rep = {}
    for sample in tqdm(samples, desc=f"loading mzML {conc}", ncols=90):
        spectra_by_rep[sample["rep"]] = load_mzml_spectra(sample, mzml_dir)

    records = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"features {conc}", ncols=90):
        if not finite(row["_aligned_rt"]) or not finite(row["_aligned_mz"]):
            continue

        fid = row.get(fid_col, idx + 1)
        aligned_rt = float(row["_aligned_rt"])
        aligned_mz = float(row["_aligned_mz"])
        original_ints = {}
        original_detected = {}

        for sample in samples:
            rep = sample["rep"]
            original_int = to_float(row.get(sample["intensity_col"]))
            original_ints[rep] = original_int
            original_detected[rep] = intensity_is_detected(original_int)

        detection_count = int(sum(original_detected.values()))
        rec = {
            "feature_ID": fid,
            "RT": aligned_rt,
            "m/z": aligned_mz,
            "detection_count": detection_count,
            "num_replicates": n_reps,
            "p_detection": assign_p_detection(detection_count, n_reps),
            "not_detected_rep_original": ";".join(str(r) for r in reps if not original_detected[r]),
        }
        if group_col:
            rec["group_ID"] = row.get(group_col)

        for rep in reps:
            rec[f"original_int_{rep}"] = original_ints[rep]
            rec[f"original_detected_{rep}"] = int(original_detected[rep])

        per_sample_windows = {}
        per_sample_csv_area = {}
        for sample in samples:
            rt_start, rt_end, window_source, csv_area = get_window_from_sample_csv(row, sample)
            per_sample_windows[sample["rep"]] = (rt_start, rt_end, window_source)
            per_sample_csv_area[sample["rep"]] = csv_area
        per_sample_windows = harmonize_windows(per_sample_windows)

        for sample in samples:
            rep = sample["rep"]
            spectra_list = spectra_by_rep.get(rep, [])
            rt_start, rt_end, window_source = per_sample_windows.get(rep, (aligned_rt - RT_TOL, aligned_rt + RT_TOL, "aligned_rt_tol"))
            rt_start, rt_end = expand_window_if_peakshape_edge(row, sample, rt_start, rt_end)
            rt_start, rt_end, window_source = clip_window_to_mzml(rt_start, rt_end, aligned_rt, spectra_list, window_source)

            all_initial, low_initial, high_initial = extract_scans_within_rt_window(spectra_list, rt_start, rt_end, aligned_mz)
            duplicate_scan_idxs, duplicate_notes = find_duplicate_scan_idxs(df, idx, sample, aligned_rt, aligned_mz, all_initial)
            all_scans, low_scans, high_scans = extract_scans_within_rt_window(
                spectra_list, rt_start, rt_end, aligned_mz, exclude_scan_idxs=duplicate_scan_idxs
            )
            target_scans = sorted(high_scans if high_scans else low_scans, key=lambda x: x["rt"])
            attrs = calculate_trace_attributes(target_scans, rt_start, rt_end)
            raw_max = attrs["raw_max"]
            raw_area = attrs["raw_area"]
            table_int = original_ints[rep]
            csv_area = per_sample_csv_area.get(rep, np.nan)

            if finite(table_int) and table_int > 0:
                resolved_int = float(table_int)
                int_source = "masscube_aligned_int"
            elif raw_max > 0:
                resolved_int = float(raw_max)
                int_source = "raw_fallback"
            else:
                resolved_int = 0.0
                int_source = "none"

            if finite(csv_area) and csv_area > 0:
                resolved_area = float(csv_area)
                area_source = "masscube_sample_csv_area"
            elif raw_area > 0:
                resolved_area = float(raw_area)
                area_source = "raw_trapezoid_fallback"
            else:
                resolved_area = 0.0
                area_source = "none"

            rec[f"int_{rep}"] = resolved_int
            rec[f"area_{rep}"] = resolved_area
            rec[f"width_{rep}"] = attrs["width"]
            rec[f"scancount_{rep}"] = attrs["scancount"]
            rec[f"smoothness_{rep}"] = attrs["smoothness"]
            rec[f"sharpness_{rep}"] = attrs["sharpness"]
            rec[f"symmetry_{rep}"] = attrs["symmetry"]
            rec[f"apex_rt_{rep}"] = attrs["apex_rt"]
            rec[f"context_int_{rep}"] = raw_max if raw_max > 0 else 0.0
            rec[f"raw_max_{rep}"] = raw_max
            rec[f"raw_area_{rep}"] = raw_area
            rec[f"rt_start_{rep}"] = rt_start
            rec[f"rt_end_{rep}"] = rt_end
            rec[f"window_source_{rep}"] = window_source
            rec[f"int_source_{rep}"] = int_source
            rec[f"area_source_{rep}"] = area_source
            rec[f"excluded_low_{rep}"] = len(low_initial)
            rec[f"excluded_dup_{rep}"] = len(duplicate_scan_idxs)

            if WRITE_SCANLIST_COLUMNS:
                rec[f"scanid_{rep}"] = scanid_string(target_scans)
                rec[f"scanlist_{rep}"] = scanlist_string(target_scans)

        records.append(rec)

    records = add_density_rank_rtshift(records, reps)
    noise_map = calculate_noise_map(records, load_blank_spectra_for_noise(BLANK_MZML_PATH))

    for rec in records:
        noise = float(noise_map.get(rec["feature_ID"], 1.0))
        rec["blank_noise_sd"] = noise
        for rep in reps:
            intensity = to_float(rec.get(f"int_{rep}"), default=0.0)
            rec[f"snr_{rep}"] = float(intensity) / noise if noise > 0 else 0.0

    full_df = safe_sort_for_feature_id(pd.DataFrame(records))
    clean_df = safe_sort_for_feature_id(build_clean_attribute_table(full_df, reps))
    analysis_df = safe_sort_for_feature_id(build_extraction_analysis_table(full_df, reps))

    clean_df.to_csv(output_dir / "masscube_pdet_full.csv", index=False)
    analysis_df.to_csv(output_dir / "masscube_extraction_analysis.csv", index=False)
    write_summary(clean_df, output_dir)
    write_attribute_files(clean_df, reps, output_dir)

    for pdet, bin_df in clean_df.groupby("p_detection", dropna=False):
        if pd.isna(pdet):
            continue
        bin_dir = output_dir / str(int(pdet))
        bin_dir.mkdir(parents=True, exist_ok=True)
        bin_df_sorted = safe_sort_for_feature_id(bin_df)
        bin_df_sorted.to_csv(bin_dir / "masscube_pdet_full.csv", index=False)
        bin_analysis = analysis_df.loc[analysis_df["feature_ID"].isin(bin_df_sorted["feature_ID"])]
        safe_sort_for_feature_id(bin_analysis).to_csv(bin_dir / "masscube_extraction_analysis.csv", index=False)

    print(f"Main table: {output_dir / 'masscube_pdet_full.csv'}")
    print(f"Diagnostics: {output_dir / 'masscube_extraction_analysis.csv'}")


def main():
    args = parse_args()
    configure_paths(args)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    if CONVERT_SINGLE_TXT_TO_CSV and not args.skip_txt_conversion:
        convert_single_txt_files_to_csv(
            SINGLE_TXT_DIR,
            SINGLE_CSV_DIR,
            overwrite=(OVERWRITE_CONVERTED_CSV or args.overwrite_converted_csv),
        )

    df, cols = prepare_feature_table(INPUT_CSV)
    sample_groups = discover_masscube_samples(df)
    if not sample_groups:
        raise ValueError("No MassCube sample columns found. Expected names like '100-1_...'.")

    print("Discovered sample groups:")
    for conc, samples in sorted(sample_groups.items()):
        print(f"  {conc}: {len(samples)} samples")

    if len(sample_groups) == 1:
        conc, samples = next(iter(sorted(sample_groups.items())))
        process_concentration(df, cols, conc, samples, OUTPUT_ROOT)
    else:
        for conc, samples in sorted(sample_groups.items()):
            process_concentration(df, cols, conc, samples, OUTPUT_ROOT / str(conc))


if __name__ == "__main__":
    main()

import argparse
from collections import defaultdict
from pathlib import Path
import re

import numpy as np
import pandas as pd
from pyteomics import mzml
from tqdm import tqdm


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

PROJECT_DIR = REPO_ROOT / "sample_data" / "mzmine_sample"
INPUT_CSV = PROJECT_DIR / "aligned_feature_table.csv"
MZML_DIR = PROJECT_DIR / "MZML" / "100"
BLANK_MZML_PATH = MZML_DIR / "MB_P1-A-4_01_13240.mzML"
OUTPUT_ROOT = PROJECT_DIR / "mzmine_pdet"

RT_COL = "rt"
MZ_COL = "mz"
FEATURE_ID_COL = "id"
GROUP_ID_COL = None
GLOBAL_RT_MIN_COL = "rt_range:min"
GLOBAL_RT_MAX_COL = "rt_range:max"

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

MZML_PATH_CACHE = {}
SPECTRA_CACHE = {}

ATTRIBUTE_NAMES = [
    "int", "area", "width", "scancount", "snr", "smoothness", "sharpness",
    "symmetry", "density", "rank", "rtshift",
]


def resolve_path(path_value, base=REPO_ROOT):
    path = Path(path_value)
    return path if path.is_absolute() else (base / path).resolve()


def configure_paths(args):
    global PROJECT_DIR, INPUT_CSV, MZML_DIR, BLANK_MZML_PATH, OUTPUT_ROOT
    global FEATURE_ID_COL, GROUP_ID_COL

    if args.project_dir:
        PROJECT_DIR = resolve_path(args.project_dir)

    INPUT_CSV = resolve_path(args.input_csv) if args.input_csv else PROJECT_DIR / "aligned_feature_table.csv"
    MZML_DIR = resolve_path(args.mzml_dir) if args.mzml_dir else PROJECT_DIR / "MZML" / "100"
    BLANK_MZML_PATH = resolve_path(args.blank_mzml) if args.blank_mzml else MZML_DIR / "MB_P1-A-4_01_13240.mzML"
    OUTPUT_ROOT = resolve_path(args.output_root) if args.output_root else PROJECT_DIR / "mzmine_pdet"

    FEATURE_ID_COL = args.feature_id_col
    GROUP_ID_COL = args.group_id_col


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate MZmine P_detection peak attributes.")
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--input-csv", default=None)
    parser.add_argument("--mzml-dir", default=None)
    parser.add_argument("--blank-mzml", default=None)
    parser.add_argument("--output-root", default=None)
    parser.add_argument("--feature-id-col", default=FEATURE_ID_COL)
    parser.add_argument("--group-id-col", default=GROUP_ID_COL)
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


def is_intensity_match(scan_int, table_int):
    if not finite(scan_int) or not finite(table_int) or float(table_int) <= 0:
        return False
    rel_diff = abs(float(scan_int) - float(table_int)) / float(table_int)
    return rel_diff <= DUP_REL_TOL or abs(float(scan_int) - float(table_int)) <= ABS_EPS


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


def discover_mzmine_samples(aligned_feature_table):
    pattern = re.compile(r"^datafile:(?P<sample>.+?\.d):(?P<field>.+)$", flags=re.IGNORECASE)
    sample_fields = defaultdict(dict)

    for col in aligned_feature_table.columns:
        match = pattern.match(str(col))
        if match:
            sample = match.group("sample")
            field = match.group("field").lower()
            sample_fields[sample][field] = col

    samples = []
    sample_pattern = re.compile(r"^(?P<conc>\d+)-(?P<rep>\d+)_")

    for sample_d, fields in sample_fields.items():
        match = sample_pattern.match(sample_d)
        if match:
            sample_base = sample_d[:-2] if sample_d.lower().endswith(".d") else sample_d
            samples.append({
                "conc": int(match.group("conc")),
                "rep": int(match.group("rep")),
                "sample_d": sample_d,
                "sample_base": sample_base,
                "height_col": fields.get("height"),
                "area_col": fields.get("area"),
                "rt_min_col": fields.get("rt_range:min"),
                "rt_max_col": fields.get("rt_range:max"),
            })

    return sorted(samples, key=lambda s: (s["conc"], s["rep"]))


def prepare_feature_table(input_csv):
    aligned_feature_table = pd.read_csv(input_csv, dtype=str, encoding="latin1")

    if RT_COL not in aligned_feature_table.columns or MZ_COL not in aligned_feature_table.columns:
        raise ValueError(f"Aligned table must contain '{RT_COL}' and '{MZ_COL}' columns.")

    aligned_feature_table["_aligned_rt"] = pd.to_numeric(aligned_feature_table[RT_COL], errors="coerce")
    aligned_feature_table["_aligned_mz"] = pd.to_numeric(aligned_feature_table[MZ_COL], errors="coerce")

    if FEATURE_ID_COL not in aligned_feature_table.columns:
        aligned_feature_table["_feature_id"] = np.arange(1, len(aligned_feature_table) + 1).astype(str)
        feature_id_col = "_feature_id"
    else:
        feature_id_col = FEATURE_ID_COL

    group_id_col = GROUP_ID_COL if GROUP_ID_COL in aligned_feature_table.columns else None
    global_rt_min_col = GLOBAL_RT_MIN_COL if GLOBAL_RT_MIN_COL in aligned_feature_table.columns else None
    global_rt_max_col = GLOBAL_RT_MAX_COL if GLOBAL_RT_MAX_COL in aligned_feature_table.columns else None

    columns = {
        "feature_id_col": feature_id_col,
        "group_id_col": group_id_col,
        "global_rt_min_col": global_rt_min_col,
        "global_rt_max_col": global_rt_max_col,
    }

    return aligned_feature_table, columns


def normalize_stem(name):
    stem = Path(name).stem
    if stem.lower().endswith(".d"):
        stem = stem[:-2]
    return stem.lower()


def find_mzml_path_for_sample(sample, mzml_dir):
    key = (sample["sample_base"].lower(), str(mzml_dir))
    if key in MZML_PATH_CACHE:
        return MZML_PATH_CACHE[key]

    mzml_dir = Path(mzml_dir)
    mzml_files = sorted(p for p in mzml_dir.iterdir() if p.is_file() and p.suffix.lower() == ".mzml")
    target = sample["sample_base"].lower()

    exact_matches = [p for p in mzml_files if normalize_stem(p.name) == target]
    if len(exact_matches) == 1:
        MZML_PATH_CACHE[key] = exact_matches[0]
        return exact_matches[0]

    prefix = f"{sample['conc']}-{sample['rep']}_".lower()
    prefix_matches = [p for p in mzml_files if p.name.lower().startswith(prefix)]
    if len(prefix_matches) == 1:
        MZML_PATH_CACHE[key] = prefix_matches[0]
        return prefix_matches[0]

    raise FileNotFoundError(f"Could not uniquely match mzML for {sample['sample_d']} in {mzml_dir}.")


def get_scan_rt_minutes(spectrum):
    scan = spectrum["scanList"]["scan"][0]
    rt = float(scan["scan start time"])
    unit = str(scan.get("unitName", "minute")).lower()
    if "second" in unit:
        rt = rt / 60.0
    return rt


def load_mzml_file(path):
    path = Path(path)
    if path in SPECTRA_CACHE:
        return SPECTRA_CACHE[path]

    spectra = []
    with mzml.read(str(path)) as reader:
        for scan_idx, spectrum in enumerate(reader):
            if int(spectrum["ms level"]) != 1:
                continue
            mzs = np.asarray(spectrum["m/z array"], dtype=np.float32)
            ints = np.asarray(spectrum["intensity array"], dtype=np.float32)
            if mzs.size == 0:
                continue
            spectra.append({
                "rt": get_scan_rt_minutes(spectrum),
                "mzs": mzs,
                "ints": ints,
                "scan_idx": scan_idx,
            })

    spectra = sorted(spectra, key=lambda s: s["rt"])
    SPECTRA_CACHE[path] = spectra
    return spectra


def load_mzml_spectra(sample, mzml_dir):
    return load_mzml_file(find_mzml_path_for_sample(sample, mzml_dir))


def get_initial_rt_window(row, sample, columns):
    aligned_rt = float(row["_aligned_rt"])

    rt_start = to_float(row.get(sample["rt_min_col"])) if sample["rt_min_col"] else np.nan
    rt_end = to_float(row.get(sample["rt_max_col"])) if sample["rt_max_col"] else np.nan

    if finite(rt_start) and finite(rt_end) and float(rt_end) > float(rt_start):
        source = "sample_rt_range"
        rt_start = float(rt_start)
        rt_end = float(rt_end)
    else:
        global_min_col = columns["global_rt_min_col"]
        global_max_col = columns["global_rt_max_col"]
        rt_start = to_float(row.get(global_min_col)) if global_min_col else np.nan
        rt_end = to_float(row.get(global_max_col)) if global_max_col else np.nan

        if finite(rt_start) and finite(rt_end) and float(rt_end) > float(rt_start):
            source = "aligned_rt_range"
            rt_start = float(rt_start)
            rt_end = float(rt_end)
        else:
            source = "aligned_rt_tol"
            rt_start = aligned_rt - RT_TOL
            rt_end = aligned_rt + RT_TOL

    if not (rt_start <= aligned_rt <= rt_end):
        rt_start = min(rt_start, aligned_rt - RT_TOL)
        rt_end = max(rt_end, aligned_rt + RT_TOL)
        source = source + "+aligned_rt_expand"

    return rt_start, rt_end, source


def harmonize_windows(sample_windows):
    good_windows = [
        (key, float(start), float(end))
        for key, (start, end, source) in sample_windows.items()
        if finite(start) and finite(end) and str(source).startswith("sample_rt_range")
    ]

    out = dict(sample_windows)

    if len(good_windows) >= 2:
        mean_start = float(np.mean([x[1] for x in good_windows]))
        mean_end = float(np.mean([x[2] for x in good_windows]))
        for key, (start, end, source) in out.items():
            if finite(start) and finite(end):
                if float(start) < mean_start - WINDOW_DIFF_THRESH or float(end) > mean_end + WINDOW_DIFF_THRESH:
                    out[key] = (mean_start, mean_end, "trimmed_to_mean")
        return out

    valid = [(float(start), float(end)) for start, end, _ in out.values() if finite(start) and finite(end)]
    if len(valid) <= 1:
        return out

    starts = np.array([x[0] for x in valid], dtype=float)
    ends = np.array([x[1] for x in valid], dtype=float)
    inconsistent = (starts.max() - starts.min() > WINDOW_DIFF_THRESH) or (ends.max() - ends.min() > WINDOW_DIFF_THRESH)

    if inconsistent:
        start = float(starts.max())
        end = float(ends.min())
        if end > start:
            for key in out:
                out[key] = (start, end, "harmonized_intersection")

    return out


def clip_window_to_mzml(rt_start, rt_end, aligned_rt, spectra, source):
    if not spectra:
        return rt_start, rt_end, source

    mzml_rt_min = spectra[0]["rt"]
    mzml_rt_max = spectra[-1]["rt"]

    if rt_end < mzml_rt_min or rt_start > mzml_rt_max:
        rt_start = max(mzml_rt_min, aligned_rt - RT_TOL)
        rt_end = min(mzml_rt_max, aligned_rt + RT_TOL)
        source = source + "+mzml_clip"

    return rt_start, rt_end, source


def extract_scans_within_rt_window(spectra, rt_start, rt_end, mz_target, excluded_scan_idxs=None):
    if not spectra or rt_start is None or rt_end is None or rt_end < rt_start:
        return [], [], []

    rts = np.array([s["rt"] for s in spectra], dtype=float)
    start_idx = np.searchsorted(rts, rt_start, side="left")
    end_idx = np.searchsorted(rts, rt_end, side="right")
    excluded_scan_idxs = excluded_scan_idxs or set()
    all_scans = []

    for spectrum in spectra[start_idx:end_idx]:
        if spectrum["scan_idx"] in excluded_scan_idxs:
            continue
        mask = mz_match_mask(spectrum["mzs"], mz_target)
        if np.any(mask):
            all_scans.append({
                "scan_idx": spectrum["scan_idx"],
                "rt": spectrum["rt"],
                "max_int": float(np.max(spectrum["ints"][mask])),
            })

    low_scans = [s for s in all_scans if not intensity_is_detected(s["max_int"])]
    high_scans = [s for s in all_scans if intensity_is_detected(s["max_int"])]
    return all_scans, low_scans, high_scans


def find_duplicate_scan_idxs(aligned_feature_table, current_idx, sample, aligned_rt, aligned_mz, scans):
    height_col = sample["height_col"]
    if height_col is None or height_col not in aligned_feature_table.columns:
        return set()

    rt_mask = aligned_feature_table["_aligned_rt"].between(aligned_rt - RT_DUP_TOL, aligned_rt + RT_DUP_TOL)
    if USE_PPM:
        ppm_diffs = np.abs((aligned_feature_table["_aligned_mz"] - aligned_mz) / aligned_mz * 1e6)
        mz_mask = ppm_diffs <= MZ_TOL_PPM
    else:
        mz_mask = aligned_feature_table["_aligned_mz"].between(aligned_mz - MZ_TOL, aligned_mz + MZ_TOL)

    nearby_features = aligned_feature_table.loc[rt_mask & mz_mask].drop(index=current_idx, errors="ignore").copy()
    if nearby_features.empty:
        return set()

    nearby_features["_candidate_height"] = pd.to_numeric(nearby_features[height_col], errors="coerce").fillna(0)
    nearby_features = nearby_features.loc[nearby_features["_candidate_height"] > 0]
    if nearby_features.empty:
        return set()

    duplicate_scan_idxs = set()
    for scan in scans:
        matched = nearby_features["_candidate_height"].apply(lambda x: is_intensity_match(scan["max_int"], x))
        if matched.any():
            duplicate_scan_idxs.add(scan["scan_idx"])

    return duplicate_scan_idxs


def scanid_string(scans):
    if not scans:
        return ""
    return ", ".join(str(s["scan_idx"] + 1) for s in sorted(scans, key=lambda x: x["rt"]))


def scanlist_string(scans):
    if not scans:
        return ""
    return "; ".join(
        f"{s['scan_idx'] + 1}|{s['rt']:.4f}|{int(round(s['max_int']))}"
        for s in sorted(scans, key=lambda x: x["rt"])
    )


def interp_rt_at_height(rts, ints, target_height):
    order = np.argsort(ints)
    sorted_ints = ints[order]
    sorted_rts = rts[order]
    unique_ints, unique_idx = np.unique(sorted_ints, return_index=True)
    unique_rts = sorted_rts[unique_idx]
    if len(unique_ints) < 2:
        return float(unique_rts[0])
    target_height = min(max(float(target_height), float(unique_ints.min())), float(unique_ints.max()))
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
    left_dist = apex_rt - rt_left
    right_dist = rt_right - apex_rt
    symmetry = right_dist / left_dist if left_dist > 0 else 1.0

    sigma = max(1.0, n_total / 6.0)
    deltas = np.diff(all_ints)
    numerator = 0.0
    denominator = 0.0

    for i, delta in enumerate(deltas):
        dist = abs(i - apex_idx)
        w_dist = np.exp(-(dist ** 2) / (2 * sigma ** 2))
        w_int = np.clip(max(all_ints[i], all_ints[i + 1]) / apex_val, 0.1, 1.0)
        abs_delta = abs(delta) * w_dist * w_int
        denominator += abs_delta
        if i < apex_idx and delta > 0:
            numerator += abs_delta
        if i >= apex_idx and delta < 0:
            numerator += abs_delta

    smoothness = numerator / denominator if denominator > 0 else 0.0
    sqrt_apex = np.sqrt(apex_val)
    sharp_vals = [
        abs(apex_val - all_ints[i]) / (abs(apex_idx - i) * sqrt_apex)
        for i in range(n_total)
        if i != apex_idx
    ]
    sharpness = float(np.max(sharp_vals)) if sharp_vals else 0.0

    return float(smoothness), float(sharpness), float(symmetry)


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


def load_blank_data(blank_path):
    if not Path(blank_path).exists():
        print(f"Blank mzML not found: {blank_path}. S/N will use noise = 1.0.")
        return []

    print(f"Reading blank mzML: {blank_path}")
    blank_data = []
    with mzml.read(str(blank_path)) as reader:
        for spectrum in reader:
            if int(spectrum["ms level"]) == 1:
                blank_data.append({
                    "mzs": np.asarray(spectrum["m/z array"], dtype=float),
                    "ints": np.asarray(spectrum["intensity array"], dtype=float),
                })

    print(f"Loaded {len(blank_data)} blank MS1 scans.")
    return blank_data


def calculate_noise_map(records, blank_data):
    noise_map = {}

    if not blank_data:
        for record in records:
            noise_map[record["feature_ID"]] = 1.0
        return noise_map

    n_nonzero = 0
    n_defaulted = 0

    for record in tqdm(records, desc="blank noise", ncols=90):
        target_mz = float(record["m/z"])
        points = []

        for spectrum in blank_data:
            mask = mz_match_mask(spectrum["mzs"], target_mz)
            points.append(float(np.max(spectrum["ints"][mask])) if np.any(mask) else 0.0)

        noise = float(np.std(points)) if points else 0.0
        if noise > 0:
            noise_map[record["feature_ID"]] = noise
            n_nonzero += 1
        else:
            noise_map[record["feature_ID"]] = 1.0
            n_defaulted += 1

    print(f"Blank noise summary: {n_nonzero} nonzero SD, {n_defaulted} defaulted to 1.0.")
    return noise_map


def add_density_rank_rtshift(records, reps):
    if not records:
        return records

    aligned_rts = np.array([float(r["RT"]) for r in records], dtype=float)

    for rep in reps:
        apex_rts = np.array([to_float(r.get(f"apex_rt_{rep}"), default=0.0) for r in records], dtype=float)
        context_ints = np.array([to_float(r.get(f"context_int_{rep}"), default=0.0) for r in records], dtype=float)

        for i, record in enumerate(records):
            curr_rt = apex_rts[i]
            curr_int = context_ints[i]

            if curr_rt <= 0 or curr_int <= 0:
                record[f"density_{rep}"] = 0
                record[f"rank_{rep}"] = 0
                record[f"rtshift_{rep}"] = 0.0
                continue

            mask = (apex_rts >= curr_rt - RT_TOL) & (apex_rts <= curr_rt + RT_TOL) & (context_ints > 0)
            local_ints = context_ints[mask]
            record[f"density_{rep}"] = int(len(local_ints))
            record[f"rank_{rep}"] = int((local_ints > curr_int).sum() + 1)
            record[f"rtshift_{rep}"] = float(curr_rt - aligned_rts[i])

    return records


def build_clean_attribute_table(full_table, reps):
    meta_cols = [c for c in ["group_ID", "feature_ID", "RT", "m/z", "p_detection"] if c in full_table.columns]
    ordered_cols = meta_cols.copy()
    for attr in ATTRIBUTE_NAMES:
        ordered_cols.extend([f"{attr}_{rep}" for rep in reps if f"{attr}_{rep}" in full_table.columns])
    return full_table[ordered_cols].copy()


def build_extraction_analysis_table(full_table, reps):
    meta_cols = [c for c in [
        "group_ID", "feature_ID", "RT", "m/z", "p_detection",
        "detection_count", "num_replicates", "not_detected_rep_original", "blank_noise_sd",
    ] if c in full_table.columns]

    prefixes = [
        "original_height", "original_detected", "int", "area", "scancount",
        "apex_rt", "context_int", "raw_max", "raw_area", "rt_start", "rt_end",
        "window_source", "int_source", "area_source", "excluded_low", "excluded_dup",
        "scanid", "scanlist",
    ]

    cols = meta_cols.copy()
    for prefix in prefixes:
        cols.extend([f"{prefix}_{rep}" for rep in reps if f"{prefix}_{rep}" in full_table.columns])
    return full_table[cols].copy()


def write_attribute_files(clean_table, reps, output_root):
    attr_dir = Path(output_root) / "attributes"
    attr_dir.mkdir(parents=True, exist_ok=True)
    meta_cols = [c for c in ["group_ID", "feature_ID", "RT", "m/z", "p_detection"] if c in clean_table.columns]

    for attr in ATTRIBUTE_NAMES:
        cols = meta_cols + [f"{attr}_{rep}" for rep in reps if f"{attr}_{rep}" in clean_table.columns]
        if len(cols) > len(meta_cols):
            clean_table[cols].to_csv(attr_dir / f"{attr}.csv", index=False)


def write_summary(clean_table, output_root):
    summary = clean_table.groupby("p_detection", dropna=False).size().reset_index(name="feature_count")
    summary = summary.sort_values("p_detection")
    summary.to_csv(Path(output_root) / "mzmine_pdet_summary.csv", index=False)


def process_aligned_feature_table(aligned_feature_table, columns):
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    strip_internal_columns(aligned_feature_table).to_csv(OUTPUT_ROOT / "original_table.csv", index=False)

    samples = discover_mzmine_samples(aligned_feature_table)
    if not samples:
        raise ValueError("No MZmine datafile:* sample columns were found.")

    reps = [s["rep"] for s in samples]
    n_reps = len(samples)

    print(f"Processing {n_reps} MZmine samples")
    print(f"mzML directory: {MZML_DIR}")
    print(f"Output directory: {OUTPUT_ROOT}")

    spectra_by_rep = {}
    for sample in tqdm(samples, desc="loading mzML", ncols=90):
        spectra_by_rep[sample["rep"]] = load_mzml_spectra(sample, MZML_DIR)

    records = []

    for idx, row in tqdm(aligned_feature_table.iterrows(), total=len(aligned_feature_table), desc="features", ncols=90):
        if not finite(row["_aligned_rt"]) or not finite(row["_aligned_mz"]):
            continue

        feature_id = row.get(columns["feature_id_col"], idx + 1)
        aligned_rt = float(row["_aligned_rt"])
        aligned_mz = float(row["_aligned_mz"])

        original_heights = {}
        original_detected = {}
        for sample in samples:
            rep = sample["rep"]
            height_col = sample["height_col"]
            original_height = to_float(row.get(height_col)) if height_col else np.nan
            original_heights[rep] = original_height
            original_detected[rep] = intensity_is_detected(original_height)

        detection_count = int(sum(original_detected.values()))
        p_detection = assign_p_detection(detection_count, n_reps)

        record = {
            "feature_ID": feature_id,
            "RT": aligned_rt,
            "m/z": aligned_mz,
            "detection_count": detection_count,
            "num_replicates": n_reps,
            "p_detection": p_detection,
            "not_detected_rep_original": ";".join([str(r) for r in reps if not original_detected[r]]),
        }

        group_col = columns["group_id_col"]
        if group_col:
            record["group_ID"] = row.get(group_col)

        for rep in reps:
            record[f"original_height_{rep}"] = original_heights[rep]
            record[f"original_detected_{rep}"] = int(original_detected[rep])

        sample_windows = {}
        for sample in samples:
            rt_start, rt_end, window_source = get_initial_rt_window(row, sample, columns)
            sample_windows[sample["rep"]] = (rt_start, rt_end, window_source)
        sample_windows = harmonize_windows(sample_windows)

        for sample in samples:
            rep = sample["rep"]
            spectra = spectra_by_rep[rep]
            rt_start, rt_end, window_source = sample_windows[rep]
            rt_start, rt_end, window_source = clip_window_to_mzml(rt_start, rt_end, aligned_rt, spectra, window_source)

            all_initial, low_initial, _ = extract_scans_within_rt_window(spectra, rt_start, rt_end, aligned_mz)
            duplicate_scan_idxs = find_duplicate_scan_idxs(aligned_feature_table, idx, sample, aligned_rt, aligned_mz, all_initial)
            _, low_scans, high_scans = extract_scans_within_rt_window(spectra, rt_start, rt_end, aligned_mz, duplicate_scan_idxs)
            target_scans = sorted(high_scans if high_scans else low_scans, key=lambda x: x["rt"])

            table_height = original_heights[rep]
            table_area = to_float(row.get(sample["area_col"])) if sample["area_col"] else np.nan
            attrs = calculate_trace_attributes(target_scans, rt_start, rt_end)

            if finite(table_height) and table_height > 0:
                resolved_int = float(table_height)
                int_source = "mzmine_height"
            elif attrs["raw_max"] > 0:
                resolved_int = attrs["raw_max"]
                int_source = "raw_fallback"
            else:
                resolved_int = 0.0
                int_source = "none"

            if finite(table_area) and table_area > 0:
                resolved_area = float(table_area)
                area_source = "mzmine_area"
            elif attrs["raw_area"] > 0:
                resolved_area = attrs["raw_area"]
                area_source = "raw_trapezoid_fallback"
            else:
                resolved_area = 0.0
                area_source = "none"

            record[f"int_{rep}"] = resolved_int
            record[f"area_{rep}"] = resolved_area
            record[f"width_{rep}"] = attrs["width"]
            record[f"scancount_{rep}"] = attrs["scancount"]
            record[f"smoothness_{rep}"] = attrs["smoothness"]
            record[f"sharpness_{rep}"] = attrs["sharpness"]
            record[f"symmetry_{rep}"] = attrs["symmetry"]
            record[f"apex_rt_{rep}"] = attrs["apex_rt"]
            record[f"context_int_{rep}"] = attrs["raw_max"] if attrs["raw_max"] > 0 else 0.0
            record[f"raw_max_{rep}"] = attrs["raw_max"]
            record[f"raw_area_{rep}"] = attrs["raw_area"]
            record[f"rt_start_{rep}"] = rt_start
            record[f"rt_end_{rep}"] = rt_end
            record[f"window_source_{rep}"] = window_source
            record[f"int_source_{rep}"] = int_source
            record[f"area_source_{rep}"] = area_source
            record[f"excluded_low_{rep}"] = len(low_initial)
            record[f"excluded_dup_{rep}"] = len(duplicate_scan_idxs)

            if WRITE_SCANLIST_COLUMNS:
                record[f"scanid_{rep}"] = scanid_string(target_scans)
                record[f"scanlist_{rep}"] = scanlist_string(target_scans)

        records.append(record)

    records = add_density_rank_rtshift(records, reps)

    blank_data = load_blank_data(BLANK_MZML_PATH)
    noise_map = calculate_noise_map(records, blank_data)
    for record in records:
        noise = noise_map[record["feature_ID"]]
        record["blank_noise_sd"] = noise
        for rep in reps:
            record[f"snr_{rep}"] = float(record.get(f"int_{rep}", 0.0)) / noise if noise > 0 else 0.0

    full_table = safe_sort_for_feature_id(pd.DataFrame(records))
    clean_table = safe_sort_for_feature_id(build_clean_attribute_table(full_table, reps))
    analysis_table = safe_sort_for_feature_id(build_extraction_analysis_table(full_table, reps))

    clean_table.to_csv(OUTPUT_ROOT / "mzmine_pdet_full.csv", index=False)
    analysis_table.to_csv(OUTPUT_ROOT / "mzmine_extraction_analysis.csv", index=False)
    write_summary(clean_table, OUTPUT_ROOT)
    write_attribute_files(clean_table, reps, OUTPUT_ROOT)

    for pdet, bin_table in clean_table.groupby("p_detection", dropna=False):
        pdet_int = int(pdet)
        bin_dir = OUTPUT_ROOT / str(pdet_int)
        bin_dir.mkdir(parents=True, exist_ok=True)
        safe_sort_for_feature_id(bin_table).to_csv(bin_dir / "mzmine_pdet_full.csv", index=False)
        bin_analysis = analysis_table.loc[analysis_table["feature_ID"].isin(bin_table["feature_ID"])]
        safe_sort_for_feature_id(bin_analysis).to_csv(bin_dir / "mzmine_extraction_analysis.csv", index=False)

    print(f"Done. Main table: {OUTPUT_ROOT / 'mzmine_pdet_full.csv'}")
    print(f"Extraction analysis: {OUTPUT_ROOT / 'mzmine_extraction_analysis.csv'}")


def main():
    args = parse_args()
    configure_paths(args)
    aligned_feature_table, columns = prepare_feature_table(INPUT_CSV)
    process_aligned_feature_table(aligned_feature_table, columns)


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mzml as pyteomics_mzml
from tqdm import tqdm


RT_TOL = 0.3
RT_DUP_TOL = 0.3
MZ_TOL = 0.01
USE_PPM = False
MZ_TOL_PPM = 20
INT_THRESHOLD = 1000
DUP_REL_TOL = 0.02
ABS_EPS = 1e-9
ALLOW_TRACE_IF_HEIGHT_MISMATCH = False
WRITE_SCANLIST_COLUMNS = True

MIN_MZ = 65.0
MAX_RT = 23.0
BLANK_INTENSITY_MULTIPLIER = 3.0
BLANK_SAMPLE_PATTERN = r"^MB"

ALIGNED_ID_COL = "Alignment ID"
ALIGNED_RT_COL = "Alignment Rt(min)"
ALIGNED_MZ_COL = "Alignment Mz"
FILL_PERCENT_COL = "Fill %"

TRACE_SCAN_COL = "Scan"
TRACE_RT_LEFT_COL = "RT left(min)"
TRACE_RT_RIGHT_COL = "RT right (min)"
TRACE_MZ_COL = "Precursor m/z"
TRACE_HEIGHT_COL = "Height"
TRACE_AREA_COL = "Area"
TRACE_SNR_COL = "S/N"

ATTRIBUTE_NAMES = [
    "int", "area", "width", "scancount", "snr", "smoothness", "sharpness",
    "symmetry", "density", "rank", "rtshift",
]


REPO_ROOT = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path.cwd()
DEFAULT_PROJECT_DIR = REPO_ROOT / "sample_data" / "msdial_sample"
DEFAULT_BLANK_INCLUDED_CSV = DEFAULT_PROJECT_DIR / "feature_table_mb.csv"
DEFAULT_NOBLANK_CSV = DEFAULT_PROJECT_DIR / "feature_table_no_mb.csv"
DEFAULT_FILTERED_OUTPUT_CSV = DEFAULT_PROJECT_DIR / "aligned_feature_table.csv"
DEFAULT_SINGLE_FILES_DIR = DEFAULT_PROJECT_DIR / "single_files"
DEFAULT_SINGLE_FILES_CSV_DIR = DEFAULT_PROJECT_DIR / "single_files_csv"
DEFAULT_MZML_DIR = DEFAULT_PROJECT_DIR / "MZML" / "100"
DEFAULT_BLANK_MZML = DEFAULT_MZML_DIR / "MB_P1-A-4_01_13240.mzML"
DEFAULT_OUTPUT_ROOT = DEFAULT_PROJECT_DIR / "msdial_pdet"

_TRACE_CACHE: dict[str, pd.DataFrame | None] = {}
_SPECTRA_CACHE: dict[Path, list[dict]] = {}
_MZML_PATH_CACHE: dict[tuple[str, Path], Path] = {}


def resolve_path(path_value: str | Path | None, repo_root: Path = REPO_ROOT) -> Path | None:
    if path_value is None:
        return None
    path = Path(path_value)
    return path if path.is_absolute() else repo_root / path


def scalar_float(x, default=np.nan):
    val = pd.to_numeric(pd.Series([x]), errors="coerce").iloc[0]
    if pd.isna(val):
        return default
    return float(val)


def finite(x) -> bool:
    val = scalar_float(x, default=np.nan)
    return bool(np.isfinite(val))


def intensity_is_detected(x) -> bool:
    return finite(x) and float(x) > INT_THRESHOLD


def mz_match_mask(mzs: np.ndarray, mz_target: float) -> np.ndarray:
    if USE_PPM:
        return np.abs((mzs - mz_target) / mz_target * 1e6) <= MZ_TOL_PPM
    return np.abs(mzs - mz_target) <= MZ_TOL


def is_intensity_match(scan_int, aligned_int) -> bool:
    if not finite(scan_int) or not finite(aligned_int) or float(aligned_int) <= 0:
        return False
    scan_int = float(scan_int)
    aligned_int = float(aligned_int)
    rel_diff = abs(scan_int - aligned_int) / aligned_int
    return rel_diff <= DUP_REL_TOL or abs(scan_int - aligned_int) <= ABS_EPS


def assign_p_detection(detection_count: int, num_reps: int) -> int:
    if num_reps <= 0:
        return 0
    return int(round(detection_count * 100.0 / float(num_reps)))


def safe_trapz(y: np.ndarray, x: np.ndarray) -> float:
    if len(y) < 2:
        return 0.0
    return float(np.trapz(y, x))


def normalize_colname(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def find_first_existing_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    norm_map = {normalize_colname(c): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        norm = normalize_colname(cand)
        if norm in norm_map:
            return norm_map[norm]
    return None


def sort_feature_table(df: pd.DataFrame) -> pd.DataFrame:
    if "feature_ID" not in df.columns:
        return df
    out = df.copy()
    out["_feature_sort"] = pd.to_numeric(out["feature_ID"], errors="coerce")
    if out["_feature_sort"].notna().any():
        out = out.sort_values(["_feature_sort", "feature_ID"], kind="mergesort")
    else:
        out = out.sort_values("feature_ID", kind="mergesort")
    return out.drop(columns=["_feature_sort"], errors="ignore")


def strip_internal_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=[c for c in df.columns if str(c).startswith("_")], errors="ignore")


def convert_trace_txts_to_csv(txt_dir: Path, csv_dir: Path, overwrite: bool = False) -> None:
    csv_dir.mkdir(parents=True, exist_ok=True)
    if not txt_dir.exists():
        print(f"No single_files directory found: {txt_dir}")
        return

    txt_files = sorted(txt_dir.glob("*.txt"))
    converted = 0
    skipped = 0

    for txt_path in tqdm(txt_files, desc="single_files txt -> csv", ncols=90):
        csv_path = csv_dir / f"{txt_path.stem}.csv"
        if csv_path.exists() and not overwrite:
            skipped += 1
            continue
        df = pd.read_csv(txt_path, sep="\t", dtype=str, encoding="latin1", engine="python")
        df.to_csv(csv_path, index=False, encoding="latin1")
        converted += 1

    print(f"Trace conversion: {converted} converted, {skipped} skipped.")


def get_alignment_columns(df: pd.DataFrame) -> dict[str, str | None]:
    return {
        "id_col": find_first_existing_col(df, [ALIGNED_ID_COL, "AlignmentID", "ID"]),
        "rt_col": find_first_existing_col(df, [ALIGNED_RT_COL, "Average Rt(min)", "Average RT(min)", "Average Rt", "RT", "Rt"]),
        "mz_col": find_first_existing_col(df, [ALIGNED_MZ_COL, "Average Mz", "Average m/z", "AverageMz", "m/z", "mz"]),
        "fill_col": find_first_existing_col(df, [FILL_PERCENT_COL, "Fill%"]),
    }


def get_sample_columns(df: pd.DataFrame) -> list[str]:
    sample_pat = re.compile(r"^\d+-\d+_")
    return [str(c) for c in df.columns if sample_pat.match(str(c))]


def get_blank_columns(df: pd.DataFrame, blank_sample_pattern: str) -> list[str]:
    pat = re.compile(blank_sample_pattern, flags=re.IGNORECASE)
    return [str(c) for c in df.columns if pat.search(str(c))]


def filter_blank_included_table(
    blank_df: pd.DataFrame,
    min_mz: float,
    max_rt: float,
    blank_multiplier: float,
    blank_sample_pattern: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = get_alignment_columns(blank_df)
    missing = [name for name, col in {"RT": cols["rt_col"], "m/z": cols["mz_col"]}.items() if col is None]
    if missing:
        raise ValueError("Missing columns in blank included table: " + ", ".join(missing))

    sample_cols = get_sample_columns(blank_df)
    blank_cols = get_blank_columns(blank_df, blank_sample_pattern)
    if not sample_cols:
        raise ValueError("No analytical sample columns were found in the blank included MS-DIAL table.")
    if not blank_cols:
        raise ValueError(f"No blank columns matched pattern {blank_sample_pattern!r}.")

    out = blank_df.copy()
    rt_values = pd.to_numeric(out[cols["rt_col"]], errors="coerce")
    mz_values = pd.to_numeric(out[cols["mz_col"]], errors="coerce")
    sample_max = out[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).max(axis=1)
    blank_max = out[blank_cols].apply(pd.to_numeric, errors="coerce").fillna(0).max(axis=1)

    # MS acquisition was set above 65 m/z, but MS-DIAL can still export lower m/z features.
    keep_mz = mz_values >= min_mz
    keep_rt = rt_values <= max_rt
    keep_blank = sample_max >= blank_multiplier * blank_max
    keep = keep_mz & keep_rt & keep_blank

    report = pd.DataFrame({
        "feature_index": np.arange(len(out), dtype=int),
        "feature_ID": out[cols["id_col"]].astype(str) if cols["id_col"] else np.arange(1, len(out) + 1).astype(str),
        "RT": rt_values,
        "m/z": mz_values,
        "max_sample_intensity": sample_max,
        "max_blank_intensity": blank_max,
        "blank_threshold": blank_multiplier * blank_max,
        "keep_mz": keep_mz,
        "keep_rt": keep_rt,
        "keep_blank": keep_blank,
        "keep": keep,
    })

    filtered = out.loc[keep].copy()
    print("Blank included table filtering")
    print(f"  input rows:                  {len(out)}")
    print(f"  retained rows:               {len(filtered)}")
    print(f"  removed m/z < {min_mz:g}:       {int((~keep_mz).sum())}")
    print(f"  removed RT > {max_rt:g} min:     {int((~keep_rt).sum())}")
    print(f"  removed sample < {blank_multiplier:g}x blank: {int((~keep_blank).sum())}")
    print(f"  analytical sample columns:    {len(sample_cols)}")
    print(f"  blank columns:                {len(blank_cols)}")

    return filtered, report


def create_filtered_noblank_table(
    noblank_csv: Path,
    blank_included_csv: Path,
    output_csv: Path,
    filtered_reference_csv: Path,
    filter_report_csv: Path,
    min_mz: float,
    max_rt: float,
    blank_multiplier: float,
    blank_sample_pattern: str,
) -> pd.DataFrame:
    blank_df = pd.read_csv(blank_included_csv, dtype=str, encoding="latin1", keep_default_na=False)
    noblank_df = pd.read_csv(noblank_csv, dtype=str, encoding="latin1", keep_default_na=False)

    reference_df, filter_report = filter_blank_included_table(
        blank_df=blank_df,
        min_mz=min_mz,
        max_rt=max_rt,
        blank_multiplier=blank_multiplier,
        blank_sample_pattern=blank_sample_pattern,
    )

    filtered_reference_csv.parent.mkdir(parents=True, exist_ok=True)
    reference_df.to_csv(filtered_reference_csv, index=False, encoding="latin1", lineterminator="\n")
    filter_report.to_csv(filter_report_csv, index=False, encoding="latin1", lineterminator="\n")

    ref_cols = get_alignment_columns(reference_df)
    new_cols = get_alignment_columns(noblank_df)
    missing = [name for name, col in {
        "blank included RT": ref_cols["rt_col"],
        "blank included m/z": ref_cols["mz_col"],
        "no blank RT": new_cols["rt_col"],
        "no blank m/z": new_cols["mz_col"],
        "no blank Fill %": new_cols["fill_col"],
    }.items() if col is None]
    if missing:
        raise ValueError("Missing columns for filtered no blank table creation: " + ", ".join(missing))

    reference_match = pd.DataFrame({
        "row_index": np.arange(len(reference_df), dtype=int),
        "rt": pd.to_numeric(reference_df[ref_cols["rt_col"]], errors="coerce"),
        "mz": pd.to_numeric(reference_df[ref_cols["mz_col"]], errors="coerce"),
    }).dropna(subset=["rt", "mz"])

    noblank_match = pd.DataFrame({
        "row_index": np.arange(len(noblank_df), dtype=int),
        "rt": pd.to_numeric(noblank_df[new_cols["rt_col"]], errors="coerce"),
        "mz": pd.to_numeric(noblank_df[new_cols["mz_col"]], errors="coerce"),
    }).dropna(subset=["rt", "mz"])

    noblank_sorted = noblank_match.sort_values(["mz", "rt", "row_index"], kind="mergesort").reset_index(drop=True)
    new_indices = noblank_sorted["row_index"].to_numpy(dtype=int)
    new_mzs = noblank_sorted["mz"].to_numpy(dtype=float)
    new_rts = noblank_sorted["rt"].to_numpy(dtype=float)

    candidates = []
    for ref in reference_match.itertuples(index=False):
        ref_index = int(ref.row_index)
        ref_rt = float(ref.rt)
        ref_mz = float(ref.mz)
        mz_window = abs(ref_mz) * MZ_TOL_PPM * 1e-6 if USE_PPM else MZ_TOL
        left = np.searchsorted(new_mzs, ref_mz - mz_window, side="left")
        right = np.searchsorted(new_mzs, ref_mz + mz_window, side="right")

        for pos in range(left, right):
            dmz = abs(float(new_mzs[pos]) - ref_mz)
            drt = abs(float(new_rts[pos]) - ref_rt)
            if drt > RT_TOL:
                continue
            if USE_PPM:
                if ref_mz == 0:
                    continue
                ppm_error = dmz / abs(ref_mz) * 1e6
                if ppm_error > MZ_TOL_PPM:
                    continue
                mz_scaled = ppm_error / MZ_TOL_PPM
            else:
                if dmz > MZ_TOL:
                    continue
                mz_scaled = dmz / MZ_TOL
            rt_scaled = drt / RT_TOL
            candidates.append({
                "reference_index": ref_index,
                "noblank_index": int(new_indices[pos]),
                "dmz": dmz,
                "drt": drt,
                "score": math.sqrt(mz_scaled ** 2 + rt_scaled ** 2),
            })

    if not candidates:
        raise ValueError("No filtered blank included features matched the no blank table.")

    candidate_df = pd.DataFrame(candidates).sort_values(
        ["score", "dmz", "drt", "reference_index", "noblank_index"],
        kind="mergesort",
    )

    used_ref = set()
    used_new = set()
    accepted = []
    for match in candidate_df.itertuples(index=False):
        ref_index = int(match.reference_index)
        new_index = int(match.noblank_index)
        if ref_index in used_ref or new_index in used_new:
            continue
        used_ref.add(ref_index)
        used_new.add(new_index)
        accepted.append({
            "reference_index": ref_index,
            "noblank_index": new_index,
            "dmz": float(match.dmz),
            "drt": float(match.drt),
            "score": float(match.score),
        })

    accepted_df = pd.DataFrame(accepted)
    filtered_df = noblank_df.iloc[sorted(accepted_df["noblank_index"].astype(int).tolist())].copy()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    filtered_df.to_csv(output_csv, index=False, encoding="latin1", lineterminator="\n")
    accepted_df.to_csv(output_csv.with_name(output_csv.stem + "_match_report.csv"), index=False, encoding="latin1", lineterminator="\n")

    print("Filtered no blank table")
    print(f"  no blank rows:               {len(noblank_df)}")
    print(f"  filtered reference rows:      {len(reference_df)}")
    print(f"  matched rows retained:        {len(filtered_df)}")
    print(f"  saved filtered table:         {output_csv}")

    return filtered_df


def prepare_feature_table(input_csv: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(input_csv, dtype=str, encoding="latin1")

    id_col = find_first_existing_col(df, [ALIGNED_ID_COL, "AlignmentID", "ID"])
    rt_col = find_first_existing_col(df, [ALIGNED_RT_COL, "Average Rt(min)", "Average RT(min)", "Average Rt", "RT", "Rt"])
    mz_col = find_first_existing_col(df, [ALIGNED_MZ_COL, "Average Mz", "Average m/z", "AverageMz", "m/z", "mz"])
    fill_col = find_first_existing_col(df, [FILL_PERCENT_COL, "Fill%"])

    missing = [name for name, col in {
        "Alignment ID": id_col,
        "RT": rt_col,
        "m/z": mz_col,
    }.items() if col is None]
    if missing:
        raise ValueError("Missing required aligned table columns: " + ", ".join(missing))

    df["_feature_id"] = df[id_col].astype(str)
    df["_aligned_rt"] = pd.to_numeric(df[rt_col], errors="coerce")
    df["_aligned_mz"] = pd.to_numeric(df[mz_col], errors="coerce")
    df["_fill_percent"] = pd.to_numeric(df[fill_col], errors="coerce") if fill_col else np.nan

    return df, {"id_col": id_col, "rt_col": rt_col, "mz_col": mz_col, "fill_col": fill_col}


def discover_msdial_samples(df: pd.DataFrame) -> dict[int, list[dict]]:
    sample_pat = re.compile(r"^(?P<conc>\d+)-(?P<rep>\d+)_")
    groups = defaultdict(list)

    for col in df.columns:
        match = sample_pat.match(str(col))
        if not match:
            continue
        conc = int(match.group("conc"))
        rep = int(match.group("rep"))
        groups[conc].append({
            "conc": conc,
            "rep": rep,
            "sample_name": str(col),
            "height_col": str(col),
        })

    return {conc: sorted(samples, key=lambda s: (s["rep"], s["sample_name"])) for conc, samples in groups.items()}


def empty_trace_df() -> pd.DataFrame:
    return pd.DataFrame(columns=[
        "_scan", "_rt_left", "_rt_right", "_prec_mz", "_height", "_area", "_snr", "_rt_center", "_width"
    ])


def standardize_trace_df(trace_df: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    col_scan = find_first_existing_col(trace_df, [TRACE_SCAN_COL])
    col_rt_left = find_first_existing_col(trace_df, [TRACE_RT_LEFT_COL, "RT left (min)", "RT left", "Rt left(min)"])
    col_rt_right = find_first_existing_col(trace_df, [TRACE_RT_RIGHT_COL, "RT right(min)", "RT right", "Rt right (min)"])
    col_mz = find_first_existing_col(trace_df, [TRACE_MZ_COL, "Precursor mz", "Precursor Mz", "m/z", "Mz"])
    col_height = find_first_existing_col(trace_df, [TRACE_HEIGHT_COL, "Peak height", "Peak Height"])
    col_area = find_first_existing_col(trace_df, [TRACE_AREA_COL, "Peak area", "Peak Area"])
    col_snr = find_first_existing_col(trace_df, [TRACE_SNR_COL, "SN", "Signal/Noise", "Signal to noise"])

    missing = [name for name, col in {
        TRACE_SCAN_COL: col_scan,
        TRACE_RT_LEFT_COL: col_rt_left,
        TRACE_RT_RIGHT_COL: col_rt_right,
        TRACE_MZ_COL: col_mz,
        TRACE_HEIGHT_COL: col_height,
        TRACE_AREA_COL: col_area,
        TRACE_SNR_COL: col_snr,
    }.items() if col is None]
    if missing:
        raise ValueError(f"Trace file {source_path} is missing columns: {', '.join(missing)}")

    df = trace_df.copy()
    df["_scan"] = pd.to_numeric(df[col_scan], errors="coerce")
    df["_rt_left"] = pd.to_numeric(df[col_rt_left], errors="coerce")
    df["_rt_right"] = pd.to_numeric(df[col_rt_right], errors="coerce")
    df["_prec_mz"] = pd.to_numeric(df[col_mz], errors="coerce")
    df["_height"] = pd.to_numeric(df[col_height], errors="coerce")
    df["_area"] = pd.to_numeric(df[col_area], errors="coerce")
    df["_snr"] = pd.to_numeric(df[col_snr], errors="coerce")
    df["_rt_center"] = (df["_rt_left"] + df["_rt_right"]) / 2.0
    df["_width"] = df["_rt_right"] - df["_rt_left"]
    return df


def trace_csv_path(sample_name: str, csv_dir: Path) -> Path | None:
    direct = csv_dir / f"{sample_name}.csv"
    if direct.exists():
        return direct
    if not csv_dir.exists():
        return None
    target = f"{sample_name}.csv".lower()
    matches = [p for p in csv_dir.iterdir() if p.name.lower() == target]
    return matches[0] if matches else None


def load_trace_csv(sample_name: str, csv_dir: Path) -> pd.DataFrame:
    cache_key = f"{csv_dir.resolve()}::{sample_name}"
    if cache_key in _TRACE_CACHE:
        cached = _TRACE_CACHE[cache_key]
        return empty_trace_df() if cached is None else cached

    path = trace_csv_path(sample_name, csv_dir)
    if path is None:
        _TRACE_CACHE[cache_key] = None
        return empty_trace_df()

    df = pd.read_csv(path, dtype=str, encoding="latin1")
    df = standardize_trace_df(df, path)
    _TRACE_CACHE[cache_key] = df
    return df


def find_matching_trace_row(trace_df: pd.DataFrame, aligned_rt: float, aligned_mz: float, aligned_height: float):
    if trace_df.empty:
        return None, "no_trace_file"

    df = trace_df.loc[trace_df["_prec_mz"].notna()].copy()
    if df.empty:
        return None, "trace_no_precursor_mz"

    if USE_PPM:
        mz_mask = np.abs((df["_prec_mz"] - aligned_mz) / aligned_mz * 1e6) <= MZ_TOL_PPM
    else:
        mz_mask = df["_prec_mz"].between(aligned_mz - MZ_TOL, aligned_mz + MZ_TOL)
    df = df.loc[mz_mask].copy()
    if df.empty:
        return None, "no_mz_match"

    inside_mask = (
        df["_rt_left"].notna()
        & df["_rt_right"].notna()
        & (df["_rt_left"] <= aligned_rt + RT_TOL)
        & (df["_rt_right"] >= aligned_rt - RT_TOL)
    )
    center_mask = df["_rt_center"].notna() & df["_rt_center"].between(aligned_rt - RT_TOL, aligned_rt + RT_TOL)
    rt_df = df.loc[inside_mask | center_mask].copy()
    if rt_df.empty:
        return None, "no_rt_match"

    rt_df["_mz_diff"] = np.abs(rt_df["_prec_mz"] - aligned_mz)
    rt_df["_rt_diff"] = np.abs(rt_df["_rt_center"].fillna(aligned_rt) - aligned_rt)

    if finite(aligned_height) and float(aligned_height) > 0:
        rt_df["_height_match"] = rt_df["_height"].apply(lambda h: is_intensity_match(h, aligned_height))
        height_matches = rt_df.loc[rt_df["_height_match"]].copy()
        if not height_matches.empty:
            height_matches = height_matches.sort_values(["_rt_diff", "_mz_diff"], kind="mergesort")
            return height_matches.iloc[0], "trace_match_rt_mz_height"
        if not ALLOW_TRACE_IF_HEIGHT_MISMATCH:
            return None, "rt_mz_match_height_mismatch"

    rt_df = rt_df.sort_values(["_rt_diff", "_mz_diff"], kind="mergesort")
    return rt_df.iloc[0], "trace_match_rt_mz_no_height_check"


def normalize_stem(path: Path) -> str:
    stem = path.stem
    if stem.lower().endswith(".d"):
        stem = stem[:-2]
    return stem.lower()


def find_mzml_path(sample: dict, mzml_dir: Path) -> Path:
    key = (sample["sample_name"].lower(), mzml_dir.resolve())
    if key in _MZML_PATH_CACHE:
        return _MZML_PATH_CACHE[key]

    files = sorted(mzml_dir.glob("*.mzML")) + sorted(mzml_dir.glob("*.mzml"))
    target = sample["sample_name"].lower()

    exact = [p for p in files if normalize_stem(p) == target]
    if len(exact) == 1:
        _MZML_PATH_CACHE[key] = exact[0]
        return exact[0]

    stem_matches = [p for p in files if p.stem.lower() in {target, target + ".d"}]
    if len(stem_matches) == 1:
        _MZML_PATH_CACHE[key] = stem_matches[0]
        return stem_matches[0]

    prefix = f"{sample['conc']}-{sample['rep']}_".lower()
    prefix_matches = [p for p in files if p.name.lower().startswith(prefix)]
    if len(prefix_matches) == 1:
        _MZML_PATH_CACHE[key] = prefix_matches[0]
        return prefix_matches[0]

    raise FileNotFoundError(f"Could not match mzML for {sample['sample_name']} in {mzml_dir}")


def get_rt_minutes(spec: dict) -> float:
    scans = spec.get("scanList", {}).get("scan", [])
    if not scans:
        return np.nan
    scan = scans[0]
    rt = scalar_float(scan.get("scan start time"), default=np.nan)
    unit = str(scan.get("scan start time unitName", "")).lower()
    if "second" in unit and np.isfinite(rt):
        rt = rt / 60.0
    return rt


def load_mzml_file(path: Path) -> list[dict]:
    path = path.resolve()
    if path in _SPECTRA_CACHE:
        return _SPECTRA_CACHE[path]

    spectra = []
    with pyteomics_mzml.read(str(path)) as reader:
        for scan_idx, spec in enumerate(reader):
            if spec.get("ms level", None) != 1:
                continue
            if "m/z array" not in spec or "intensity array" not in spec:
                continue
            rt = get_rt_minutes(spec)
            if not np.isfinite(rt):
                continue
            mzs = np.asarray(spec["m/z array"], dtype=float)
            ints = np.asarray(spec["intensity array"], dtype=float)
            if mzs.size == 0 or ints.size == 0:
                continue
            spectra.append({"rt": float(rt), "mzs": mzs, "ints": ints, "scan_idx": int(scan_idx)})

    spectra = sorted(spectra, key=lambda s: s["rt"])
    _SPECTRA_CACHE[path] = spectra
    return spectra


def load_mzml_spectra(sample: dict, mzml_dir: Path) -> list[dict]:
    return load_mzml_file(find_mzml_path(sample, mzml_dir))


def clip_window_to_mzml(rt_start: float, rt_end: float, aligned_rt: float, spectra_list: list[dict], source: str):
    if not spectra_list:
        return rt_start, rt_end, source
    mzml_rt_min = spectra_list[0]["rt"]
    mzml_rt_max = spectra_list[-1]["rt"]
    if rt_end < mzml_rt_min or rt_start > mzml_rt_max:
        rt_start = max(mzml_rt_min, aligned_rt - RT_TOL)
        rt_end = min(mzml_rt_max, aligned_rt + RT_TOL)
        source = source + "+mzml_clip"
    return rt_start, rt_end, source


def extract_scans(spectra_list: list[dict], rt_start: float, rt_end: float, mz_target: float, exclude_scan_idxs=None):
    if not spectra_list or rt_end < rt_start:
        return [], [], []

    rts = np.array([s["rt"] for s in spectra_list], dtype=float)
    start_idx = np.searchsorted(rts, rt_start, side="left")
    end_idx = np.searchsorted(rts, rt_end, side="right")
    exclude_scan_idxs = exclude_scan_idxs or set()
    all_scans = []

    for spec in spectra_list[start_idx:end_idx]:
        if spec["scan_idx"] in exclude_scan_idxs:
            continue
        mask = mz_match_mask(spec["mzs"], mz_target)
        if not np.any(mask):
            continue
        all_scans.append({
            "scan_idx": spec["scan_idx"],
            "rt": spec["rt"],
            "max_int": float(np.max(spec["ints"][mask])),
        })

    lowint = [s for s in all_scans if not intensity_is_detected(s["max_int"])]
    highint = [s for s in all_scans if intensity_is_detected(s["max_int"])]
    return all_scans, lowint, highint


def find_duplicate_scan_idxs(df, current_idx, sample, aligned_rt, aligned_mz, scans_to_check):
    height_col = sample["height_col"]
    rt_mask = df["_aligned_rt"].between(aligned_rt - RT_DUP_TOL, aligned_rt + RT_DUP_TOL)
    mz_mask = df["_aligned_mz"].between(aligned_mz - MZ_TOL, aligned_mz + MZ_TOL)
    candidates = df.loc[rt_mask & mz_mask].copy()
    if current_idx in candidates.index:
        candidates = candidates.drop(index=current_idx)
    if candidates.empty:
        return set()

    candidates["_candidate_height"] = pd.to_numeric(candidates[height_col], errors="coerce").fillna(0)
    candidates = candidates.loc[candidates["_candidate_height"] > 0]
    duplicate_scan_idxs = set()

    for scan in scans_to_check:
        for _, cand in candidates.iterrows():
            if is_intensity_match(scan["max_int"], float(cand["_candidate_height"])):
                duplicate_scan_idxs.add(scan["scan_idx"])
                break

    return duplicate_scan_idxs


def scanlist_string(scans: list[dict]) -> str:
    scans = sorted(scans, key=lambda x: x["rt"])
    return "; ".join(f"{s['scan_idx'] + 1}|{s['rt']:.4f}|{int(round(s['max_int']))}" for s in scans)


def scanid_string(scans: list[dict]) -> str:
    scans = sorted(scans, key=lambda x: x["rt"])
    return ", ".join(str(s["scan_idx"] + 1) for s in scans)


def calculate_shape_metrics(scans: list[dict]):
    if not scans:
        return 0.0, 0.0, 1.0

    scans = sorted(scans, key=lambda x: x["rt"])
    rts = np.array([float(s["rt"]) for s in scans], dtype=float)
    ints = np.array([float(s["max_int"]) for s in scans], dtype=float)
    n_total = len(ints)

    if n_total < 3:
        return 0.0, 0.0, 1.0

    apex_idx = int(np.argmax(ints))
    apex_val = float(ints[apex_idx])
    apex_rt = float(rts[apex_idx])
    if apex_val <= 0:
        return 0.0, 0.0, 1.0

    target_height = 0.10 * apex_val
    rt_left = np.interp(target_height, ints[:apex_idx + 1], rts[:apex_idx + 1])
    rt_right = np.interp(target_height, ints[apex_idx:][::-1], rts[apex_idx:][::-1])
    A = apex_rt - rt_left
    B = rt_right - apex_rt
    symmetry = B / A if A > 0 else 1.0

    sigma = max(1.0, n_total / 6.0)
    deltas = np.diff(ints)
    numerator = 0.0
    denominator = 0.0

    for i, d in enumerate(deltas):
        dist = abs(i - apex_idx)
        w_dist = np.exp(-(dist ** 2) / (2 * sigma ** 2))
        w_int = np.clip(max(ints[i], ints[i + 1]) / apex_val, 0.1, 1.0)
        weighted_delta = abs(d) * w_dist * w_int
        denominator += weighted_delta
        if i < apex_idx and d > 0:
            numerator += weighted_delta
        if i >= apex_idx and d < 0:
            numerator += weighted_delta

    smoothness = numerator / denominator if denominator > 0 else 0.0
    sharp_vals = [
        abs(apex_val - ints[i]) / (abs(apex_idx - i) * np.sqrt(apex_val))
        for i in range(n_total)
        if i != apex_idx
    ]
    sharpness = float(np.max(sharp_vals)) if sharp_vals else 0.0
    return float(smoothness), float(sharpness), float(symmetry)


def calculate_trace_attributes(scans: list[dict], rt_start: float, rt_end: float) -> dict:
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
    smoothness, sharpness, symmetry = calculate_shape_metrics(scans)

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


def load_blank_spectra(blank_path: Path) -> list[dict]:
    if not blank_path.exists():
        print(f"Blank file not found. Fallback S/N will use noise = 1: {blank_path}")
        return []
    return load_mzml_file(blank_path)


def calculate_noise_map(records: list[dict], blank_data: list[dict]) -> dict:
    if not blank_data:
        return {rec["feature_ID"]: 1.0 for rec in records}

    noise_map = {}
    nonzero = 0
    defaulted = 0

    for rec in tqdm(records, desc="blank noise", ncols=90):
        target_mz = float(rec["m/z"])
        points = []
        for spec in blank_data:
            mask = mz_match_mask(spec["mzs"], target_mz)
            points.append(float(np.max(spec["ints"][mask])) if np.any(mask) else 0.0)
        noise = float(np.std(points)) if points else 0.0
        if noise > 0:
            noise_map[rec["feature_ID"]] = noise
            nonzero += 1
        else:
            noise_map[rec["feature_ID"]] = 1.0
            defaulted += 1

    print(f"Blank noise: {nonzero} nonzero SD; {defaulted} defaulted to 1.0.")
    return noise_map


def add_density_rank_rtshift(records: list[dict], reps: list[int]) -> list[dict]:
    if not records:
        return records

    aligned_rts = np.array([float(r["RT"]) for r in records], dtype=float)
    for rep in reps:
        apex_rts = np.array([scalar_float(r.get(f"apex_rt_{rep}"), default=0.0) for r in records], dtype=float)
        context_ints = np.array([scalar_float(r.get(f"context_int_{rep}"), default=0.0) for r in records], dtype=float)

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


def pdet_from_fill_percent(row: pd.Series, n_reps: int) -> tuple[float, int, int]:
    fill_fraction = scalar_float(row.get("_fill_percent"), default=0.0)
    if fill_fraction > 1:
        fill_fraction /= 100.0
    fill_fraction = float(np.clip(fill_fraction, 0.0, 1.0))
    detection_count = int(round(fill_fraction * n_reps))
    return fill_fraction, detection_count, assign_p_detection(detection_count, n_reps)


def write_summary(clean_df: pd.DataFrame, output_dir: Path) -> None:
    summary = clean_df.groupby("p_detection", dropna=False).size().reset_index(name="feature_count")
    summary = summary.sort_values("p_detection")
    summary.to_csv(output_dir / "msdial_pdet_summary.csv", index=False)


def build_clean_table(full_df: pd.DataFrame, reps: list[int]) -> pd.DataFrame:
    meta_cols = [c for c in ["feature_ID", "RT", "m/z", "fill_percent", "p_detection"] if c in full_df.columns]
    cols = meta_cols.copy()
    for attr in ATTRIBUTE_NAMES:
        cols.extend([f"{attr}_{r}" for r in reps if f"{attr}_{r}" in full_df.columns])
    return full_df[cols].copy()


def build_analysis_table(full_df: pd.DataFrame, reps: list[int]) -> pd.DataFrame:
    meta_cols = [c for c in [
        "feature_ID", "RT", "m/z", "fill_percent", "p_detection", "detection_count",
        "num_replicates", "not_detected_rep_original", "blank_noise_sd",
    ] if c in full_df.columns]
    prefixes = [
        "original_height", "original_detected", "trace_match_status", "trace_scan", "trace_height",
        "trace_area", "trace_snr", "int", "area", "snr", "scancount", "apex_rt", "context_int",
        "raw_max", "raw_area", "rt_start", "rt_end", "window_source", "int_source", "area_source",
        "snr_source", "excluded_low", "excluded_dup", "scanid", "scanlist",
    ]
    cols = meta_cols.copy()
    for prefix in prefixes:
        cols.extend([f"{prefix}_{r}" for r in reps if f"{prefix}_{r}" in full_df.columns])
    return full_df[[c for c in cols if c in full_df.columns]].copy()


def write_attribute_files(clean_df: pd.DataFrame, reps: list[int], output_dir: Path) -> None:
    attr_dir = output_dir / "attributes"
    attr_dir.mkdir(parents=True, exist_ok=True)
    meta_cols = [c for c in ["feature_ID", "RT", "m/z", "fill_percent", "p_detection"] if c in clean_df.columns]
    for attr in ATTRIBUTE_NAMES:
        cols = meta_cols + [f"{attr}_{r}" for r in reps if f"{attr}_{r}" in clean_df.columns]
        if len(cols) > len(meta_cols):
            clean_df[cols].to_csv(attr_dir / f"{attr}.csv", index=False)


def process_concentration(
    df: pd.DataFrame,
    conc: int,
    samples: list[dict],
    mzml_dir: Path,
    trace_csv_dir: Path,
    blank_mzml: Path,
    output_dir: Path,
    use_fill_percent_pdet: bool,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    strip_internal_columns(df).to_csv(output_dir / "original_table.csv", index=False)

    reps = [s["rep"] for s in samples]
    n_reps = len(samples)

    print(f"\nProcessing MS-DIAL concentration {conc}")
    print(f"mzML: {mzml_dir}")
    print(f"output: {output_dir}")

    spectra_by_rep = {sample["rep"]: load_mzml_spectra(sample, mzml_dir) for sample in tqdm(samples, desc="loading mzML", ncols=90)}
    records = []

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="features", ncols=90):
        if not finite(row["_aligned_rt"]) or not finite(row["_aligned_mz"]):
            continue

        fid = row["_feature_id"]
        aligned_rt = float(row["_aligned_rt"])
        aligned_mz = float(row["_aligned_mz"])

        original_heights = {sample["rep"]: scalar_float(row.get(sample["height_col"]), default=np.nan) for sample in samples}
        original_detected = {rep: intensity_is_detected(val) for rep, val in original_heights.items()}

        fill_fraction = scalar_float(row.get("_fill_percent"), default=np.nan)
        if use_fill_percent_pdet and np.isfinite(fill_fraction):
            fill_fraction, detection_count, p_detection = pdet_from_fill_percent(row, n_reps)
        else:
            detection_count = int(sum(original_detected.values()))
            p_detection = assign_p_detection(detection_count, n_reps)
            if np.isfinite(fill_fraction) and fill_fraction > 1:
                fill_fraction = fill_fraction / 100.0

        rec = {
            "feature_ID": fid,
            "RT": aligned_rt,
            "m/z": aligned_mz,
            "fill_percent": fill_fraction,
            "detection_count": detection_count,
            "num_replicates": n_reps,
            "p_detection": p_detection,
            "not_detected_rep_original": ";".join([str(r) for r in reps if not original_detected[r]]),
        }

        for rep in reps:
            rec[f"original_height_{rep}"] = original_heights[rep]
            rec[f"original_detected_{rep}"] = int(original_detected[rep])

        for sample in samples:
            rep = sample["rep"]
            aligned_height = original_heights[rep]
            trace_df = load_trace_csv(sample["sample_name"], trace_csv_dir)
            trace_row, trace_status = find_matching_trace_row(trace_df, aligned_rt, aligned_mz, aligned_height)

            trace_height = scalar_float(trace_row.get("_height"), default=np.nan) if trace_row is not None else np.nan
            trace_area = scalar_float(trace_row.get("_area"), default=np.nan) if trace_row is not None else np.nan
            trace_snr = scalar_float(trace_row.get("_snr"), default=np.nan) if trace_row is not None else np.nan
            trace_scan = scalar_float(trace_row.get("_scan"), default=np.nan) if trace_row is not None else np.nan

            if trace_row is not None and finite(trace_row.get("_rt_left")) and finite(trace_row.get("_rt_right")):
                rt_start = float(trace_row["_rt_left"])
                rt_end = float(trace_row["_rt_right"])
                window_source = "msdial_trace"
            else:
                rt_start = aligned_rt - RT_TOL
                rt_end = aligned_rt + RT_TOL
                window_source = "aligned_rt_tol"

            if not (rt_start <= aligned_rt <= rt_end):
                rt_start = min(rt_start, aligned_rt - RT_TOL)
                rt_end = max(rt_end, aligned_rt + RT_TOL)
                window_source = window_source + "+aligned_rt_expand"

            spectra_list = spectra_by_rep[rep]
            rt_start, rt_end, window_source = clip_window_to_mzml(rt_start, rt_end, aligned_rt, spectra_list, window_source)

            all_initial, low_initial, _ = extract_scans(spectra_list, rt_start, rt_end, aligned_mz)
            duplicate_scan_idxs = find_duplicate_scan_idxs(df, idx, sample, aligned_rt, aligned_mz, all_initial)
            _, low_scans, high_scans = extract_scans(spectra_list, rt_start, rt_end, aligned_mz, duplicate_scan_idxs)
            target_scans = sorted(high_scans if high_scans else low_scans, key=lambda x: x["rt"])
            attrs = calculate_trace_attributes(target_scans, rt_start, rt_end)

            raw_max = attrs["raw_max"]
            raw_area = attrs["raw_area"]

            if finite(trace_height) and trace_height > 0:
                resolved_int = float(trace_height)
                int_source = "msdial_trace_height"
            elif raw_max > 0:
                resolved_int = float(raw_max)
                int_source = "raw_fallback"
            elif finite(aligned_height) and aligned_height > 0:
                resolved_int = float(aligned_height)
                int_source = "aligned_table_fallback"
            else:
                resolved_int = 0.0
                int_source = "none"

            if finite(trace_area) and trace_area > 0:
                resolved_area = float(trace_area)
                area_source = "msdial_trace_area"
            elif raw_area > 0:
                resolved_area = float(raw_area)
                area_source = "raw_trapezoid_fallback"
            else:
                resolved_area = 0.0
                area_source = "none"

            if finite(trace_snr) and trace_snr > 0:
                resolved_snr = float(trace_snr)
                snr_source = "msdial_trace_snr"
            else:
                resolved_snr = np.nan
                snr_source = "pending_blank_noise_fallback"

            context_int = raw_max if raw_max > 0 else resolved_int

            rec[f"int_{rep}"] = resolved_int
            rec[f"area_{rep}"] = resolved_area
            rec[f"width_{rep}"] = attrs["width"]
            rec[f"scancount_{rep}"] = attrs["scancount"]
            rec[f"snr_{rep}"] = resolved_snr
            rec[f"smoothness_{rep}"] = attrs["smoothness"]
            rec[f"sharpness_{rep}"] = attrs["sharpness"]
            rec[f"symmetry_{rep}"] = attrs["symmetry"]
            rec[f"apex_rt_{rep}"] = attrs["apex_rt"]
            rec[f"context_int_{rep}"] = context_int
            rec[f"raw_max_{rep}"] = raw_max
            rec[f"raw_area_{rep}"] = raw_area
            rec[f"rt_start_{rep}"] = rt_start
            rec[f"rt_end_{rep}"] = rt_end
            rec[f"window_source_{rep}"] = window_source
            rec[f"int_source_{rep}"] = int_source
            rec[f"area_source_{rep}"] = area_source
            rec[f"snr_source_{rep}"] = snr_source
            rec[f"trace_match_status_{rep}"] = trace_status
            rec[f"trace_scan_{rep}"] = trace_scan
            rec[f"trace_height_{rep}"] = trace_height
            rec[f"trace_area_{rep}"] = trace_area
            rec[f"trace_snr_{rep}"] = trace_snr
            rec[f"excluded_low_{rep}"] = len(low_initial)
            rec[f"excluded_dup_{rep}"] = len(duplicate_scan_idxs)

            if WRITE_SCANLIST_COLUMNS:
                rec[f"scanid_{rep}"] = scanid_string(target_scans)
                rec[f"scanlist_{rep}"] = scanlist_string(target_scans)

        records.append(rec)

    records = add_density_rank_rtshift(records, reps)
    noise_map = calculate_noise_map(records, load_blank_spectra(blank_mzml))

    for rec in records:
        noise = noise_map.get(rec["feature_ID"], 1.0)
        rec["blank_noise_sd"] = noise
        for rep in reps:
            if not finite(rec.get(f"snr_{rep}")) or rec.get(f"snr_{rep}") <= 0:
                intensity = scalar_float(rec.get(f"int_{rep}"), default=0.0)
                rec[f"snr_{rep}"] = intensity / noise if noise > 0 else 0.0
                rec[f"snr_source_{rep}"] = "blank_noise_fallback"

    full_df = sort_feature_table(pd.DataFrame(records))
    clean_df = sort_feature_table(build_clean_table(full_df, reps))
    analysis_df = sort_feature_table(build_analysis_table(full_df, reps))

    clean_df.to_csv(output_dir / "msdial_pdet_full.csv", index=False)
    analysis_df.to_csv(output_dir / "msdial_extraction_analysis.csv", index=False)
    write_summary(clean_df, output_dir)
    write_attribute_files(clean_df, reps, output_dir)

    for pdet, bin_df in clean_df.groupby("p_detection", dropna=False):
        pdet_int = int(pdet)
        bin_dir = output_dir / str(pdet_int)
        bin_dir.mkdir(parents=True, exist_ok=True)
        sort_feature_table(bin_df).to_csv(bin_dir / "msdial_pdet_full.csv", index=False)
        bin_analysis = analysis_df.loc[analysis_df["feature_ID"].isin(bin_df["feature_ID"])]
        sort_feature_table(bin_analysis).to_csv(bin_dir / "msdial_extraction_analysis.csv", index=False)

    print(f"Saved: {output_dir / 'msdial_pdet_full.csv'}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", default=str(DEFAULT_PROJECT_DIR))
    parser.add_argument("--blank-included-csv", "--blank-filtered-reference-csv", dest="blank_included_csv", default=None)
    parser.add_argument("--noblank-csv", "--noblank-unfiltered-csv", dest="noblank_csv", default=None)
    parser.add_argument("--filtered-output-csv", default=None)
    parser.add_argument("--single-files-dir", default=None)
    parser.add_argument("--single-files-csv-dir", default=None)
    parser.add_argument("--mzml-dir", default=None)
    parser.add_argument("--blank-mzml", default=None)
    parser.add_argument("--output-root", default=None)
    parser.add_argument("--overwrite-converted-csv", action="store_true")
    parser.add_argument("--skip-txt-conversion", action="store_true")
    parser.add_argument("--skip-blank-filtering", action="store_true")
    parser.add_argument("--pdet-from-intensities", action="store_true")
    parser.add_argument("--min-mz", type=float, default=MIN_MZ)
    parser.add_argument("--max-rt", type=float, default=MAX_RT)
    parser.add_argument("--blank-multiplier", type=float, default=BLANK_INTENSITY_MULTIPLIER)
    parser.add_argument("--blank-sample-pattern", default=BLANK_SAMPLE_PATTERN)
    return parser.parse_args()


def main():
    args = parse_args()
    project_dir = resolve_path(args.project_dir)
    blank_included_csv = resolve_path(args.blank_included_csv) if args.blank_included_csv else project_dir / "feature_table_mb.csv"
    noblank_csv = resolve_path(args.noblank_csv) if args.noblank_csv else project_dir / "feature_table_no_mb.csv"
    filtered_output_csv = resolve_path(args.filtered_output_csv) if args.filtered_output_csv else project_dir / "aligned_feature_table.csv"
    single_files_dir = resolve_path(args.single_files_dir) if args.single_files_dir else project_dir / "single_files"
    single_files_csv_dir = resolve_path(args.single_files_csv_dir) if args.single_files_csv_dir else project_dir / "single_files_csv"
    mzml_dir = resolve_path(args.mzml_dir) if args.mzml_dir else project_dir / "MZML" / "100"
    blank_mzml = resolve_path(args.blank_mzml) if args.blank_mzml else mzml_dir / "MB_P1-A-4_01_13240.mzML"
    output_root = resolve_path(args.output_root) if args.output_root else project_dir / "msdial_pdet"
    output_root.mkdir(parents=True, exist_ok=True)

    if args.skip_blank_filtering:
        input_csv = filtered_output_csv
    else:
        input_csv = filtered_output_csv
        create_filtered_noblank_table(
            noblank_csv=noblank_csv,
            blank_included_csv=blank_included_csv,
            output_csv=filtered_output_csv,
            filtered_reference_csv=output_root / "msdial_blank_filtered_reference.csv",
            filter_report_csv=output_root / "msdial_blank_filter_report.csv",
            min_mz=args.min_mz,
            max_rt=args.max_rt,
            blank_multiplier=args.blank_multiplier,
            blank_sample_pattern=args.blank_sample_pattern,
        )

    if not args.skip_txt_conversion:
        convert_trace_txts_to_csv(single_files_dir, single_files_csv_dir, overwrite=args.overwrite_converted_csv)

    df, _ = prepare_feature_table(input_csv)
    sample_groups = discover_msdial_samples(df)

    if not sample_groups:
        raise ValueError("No MS-DIAL sample columns were found, e.g. 100-13_P1-C-1_01_13570.")

    print("Discovered sample groups:")
    for conc, samples in sample_groups.items():
        print(f"  {conc}: {len(samples)} samples")

    use_fill_percent_pdet = not args.pdet_from_intensities
    if len(sample_groups) == 1:
        conc, samples = next(iter(sorted(sample_groups.items())))
        process_concentration(
            df, conc, samples, mzml_dir, single_files_csv_dir, blank_mzml, output_root,
            use_fill_percent_pdet=use_fill_percent_pdet,
        )
    else:
        for conc, samples in sorted(sample_groups.items()):
            process_concentration(
                df, conc, samples, mzml_dir, single_files_csv_dir, blank_mzml, output_root / str(conc),
                use_fill_percent_pdet=use_fill_percent_pdet,
            )


if __name__ == "__main__":
    main()
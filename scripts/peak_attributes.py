import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

DEFAULT_MODEL = "masscube"

MODEL_CONFIG = {
    "masscube": {
        "project_dir": REPO_ROOT / "sample_data" / "masscube_sample",
        "input_name": "masscube_pdet_full.csv",
        "pdet_folder": "masscube_pdet",
        "summary_prefix": "masscube",
    },
    "mzmine": {
        "project_dir": REPO_ROOT / "sample_data" / "mzmine_sample",
        "input_name": "mzmine_pdet_full.csv",
        "pdet_folder": "mzmine_pdet",
        "summary_prefix": "mzmine",
    },
    "msdial": {
        "project_dir": REPO_ROOT / "sample_data" / "msdial_sample",
        "input_name": "msdial_pdet_full.csv",
        "pdet_folder": "msdial_pdet",
        "summary_prefix": "msdial",
    },
}

METRICS = [
    "int", "area", "width", "scancount", "snr", "smoothness", "sharpness",
    "symmetry", "density", "rank", "rtshift",
]


def resolve_path(path_value, base=REPO_ROOT):
    path = Path(path_value)
    return path if path.is_absolute() else (base / path).resolve()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize replicate-level peak attributes from *_pdet_full.csv."
    )
    parser.add_argument(
        "--model",
        choices=sorted(MODEL_CONFIG),
        default=DEFAULT_MODEL,
        help="Default is masscube. Change to mzmine or msdial to process those outputs.",
    )
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--input-csv", default=None)
    parser.add_argument("--output-dir", default=None)
    return parser.parse_args()


def default_input_csv(model, project_dir=None):
    config = MODEL_CONFIG[model]
    root = resolve_path(project_dir) if project_dir else config["project_dir"]
    return root / config["pdet_folder"] / config["input_name"]


def infer_model_from_input(input_csv, fallback_model):
    name = input_csv.name.lower()
    parent_parts = [p.lower() for p in input_csv.parts]

    for model in MODEL_CONFIG:
        if name.startswith(model) or any(model in part for part in parent_parts):
            return model

    return fallback_model


def default_output_dir(input_csv):
    return input_csv.parent / "peak_analysis"


def metric_columns(df, metric):
    return [c for c in df.columns if c.startswith(f"{metric}_")]


def numeric_values(row, cols):
    if not cols:
        return np.array([], dtype=float)
    values = pd.to_numeric(row[cols], errors="coerce").to_numpy(dtype=float)
    return values[~np.isnan(values)]


def compute_stats(values):
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if values.size == 0:
        return {
            "values_str": "",
            "max": np.nan,
            "min": np.nan,
            "mean": np.nan,
            "median": np.nan,
            "75_perc": np.nan,
            "25_perc": np.nan,
            "sd": np.nan,
            "rsd": np.nan,
        }

    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1)) if values.size > 1 else np.nan
    rsd = float(sd / mean * 100) if values.size > 1 and mean != 0 and not np.isnan(sd) else np.nan

    return {
        "values_str": ";".join(str(v) for v in values),
        "max": float(np.max(values)),
        "min": float(np.min(values)),
        "mean": mean,
        "median": float(np.median(values)),
        "75_perc": float(np.percentile(values, 75)),
        "25_perc": float(np.percentile(values, 25)),
        "sd": sd,
        "rsd": rsd,
    }


def required_metadata_columns(df):
    required = ["feature_ID", "p_detection", "m/z", "RT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input table is missing required columns: {missing}")
    return required


def summarize_metric(df, metric, cols):
    rows = []

    for _, row in df.iterrows():
        values = numeric_values(row, cols)
        stats = compute_stats(values)

        rows.append({
            "feature_ID": row["feature_ID"],
            "p_detection": int(row["p_detection"]) if not pd.isna(row["p_detection"]) else np.nan,
            "m/z": row["m/z"],
            "RT": row["RT"],
            f"{metric}_values": stats["values_str"],
            "max": stats["max"],
            "min": stats["min"],
            "mean": stats["mean"],
            "median": stats["median"],
            "75_perc": stats["75_perc"],
            "25_perc": stats["25_perc"],
            "sd": stats["sd"],
            "rsd": stats["rsd"],
        })

    return pd.DataFrame(rows)


def write_metric_summaries(df, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    required_metadata_columns(df)

    written = []
    skipped = []

    for metric in METRICS:
        cols = metric_columns(df, metric)
        if not cols:
            skipped.append(metric)
            continue

        summary_df = summarize_metric(df, metric, cols)
        output_path = output_dir / f"{metric}_summary.csv"
        summary_df.to_csv(output_path, index=False)
        written.append(output_path)

    return written, skipped


def main():
    args = parse_args()

    input_csv = resolve_path(args.input_csv) if args.input_csv else default_input_csv(args.model, args.project_dir)
    model = infer_model_from_input(input_csv, args.model)
    output_dir = resolve_path(args.output_dir) if args.output_dir else default_output_dir(input_csv)

    if not input_csv.exists():
        raise FileNotFoundError(f"Input file not found: {input_csv}")

    print(f"Model: {model}")
    print(f"Input: {input_csv}")
    print(f"Output: {output_dir}")

    df = pd.read_csv(input_csv)
    written, skipped = write_metric_summaries(df, output_dir)

    print(f"Wrote {len(written)} summary files.")
    if skipped:
        print("Skipped metrics with no columns: " + ", ".join(skipped))
    print("Peak attribute summary complete.")


if __name__ == "__main__":
    main()
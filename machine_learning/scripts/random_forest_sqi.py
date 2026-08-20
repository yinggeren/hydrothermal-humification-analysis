"""Random Forest analysis for hydrochar properties and soil quality index.

This script reproduces the feature-importance calculation used for Fig. 5c.
The final model keeps pH and excludes application rate (AR, %).

Inputs:
    data/tableS7_hydrochar_properties.csv
    data/tableS8_soil_quality_indicators.csv

Outputs:
    outputs/sqi_dataset_for_random_forest.csv
    outputs/fig5c_feature_importance.csv
    outputs/fig5c_feature_importance.svg
    outputs/run_metadata.csv
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd


RF_FEATURES = ["pH", "C (%)", "N (%)", "P (g/kg)", "K (g/kg)", "O/C", "H/C"]
TARGET = "SQI"
SOIL_INDICATORS = [
    "pH",
    "EC (us/cm)",
    "CEC (cmol/kg)",
    "OM (g/kg)",
    "TC (mg/g)",
    "MBC mg/g",
    "TN (mg/g)",
    "TP (mg/g)",
    "TK (mg/g)",
    "NH4+-N (mg/g)",
    "NO3--N (mg/g)",
    "AN (mg/g)",
    "AP (mg/g)",
    "AK (mg/g)",
]


def parse_number(value: object) -> float:
    if value is None:
        return np.nan
    text = str(value).strip()
    if not text:
        return np.nan
    text = (
        text.replace("\u2212", "-")
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace(",", "")
    )
    if text.lower() in {"na", "n/a", "nd", "-", "--", "none"}:
        return np.nan
    import re

    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    if not match:
        return np.nan
    return float(match.group(0))


def mean_impute(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    matrix = matrix.astype(float, copy=True)
    means = np.nanmean(matrix, axis=0)
    means = np.where(np.isnan(means), 0.0, means)
    row_idx, col_idx = np.where(np.isnan(matrix))
    matrix[row_idx, col_idx] = means[col_idx]
    return matrix, means


def minmax_scale(matrix: np.ndarray) -> np.ndarray:
    mins = np.min(matrix, axis=0)
    maxs = np.max(matrix, axis=0)
    denom = np.where(maxs > mins, maxs - mins, 1.0)
    return (matrix - mins) / denom


def calculate_pca_sqi(
    soil_df: pd.DataFrame, n_components: int = 3
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    missing = [col for col in SOIL_INDICATORS if col not in soil_df.columns]
    if missing:
        raise ValueError(f"Soil table is missing required columns: {missing}")

    soil_numeric = soil_df[SOIL_INDICATORS].map(parse_number).to_numpy(dtype=float)
    soil_numeric, _ = mean_impute(soil_numeric)
    soil_scaled = minmax_scale(soil_numeric)

    centered = soil_scaled - soil_scaled.mean(axis=0)
    covariance = np.cov(centered, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    total = eigenvalues.sum()
    explained = eigenvalues / total if total > 0 else np.zeros_like(eigenvalues)
    n_components = min(n_components, len(SOIL_INDICATORS))
    retained_loadings = eigenvectors[:, :n_components]
    retained_variance = explained[:n_components]

    raw_weights = np.sum(np.abs(retained_loadings) * retained_variance, axis=1)
    weights = raw_weights / raw_weights.sum()
    sqi = soil_scaled @ weights

    weights_df = pd.DataFrame(
        {
            "indicator": SOIL_INDICATORS,
            "weight": weights,
            "raw_weight": raw_weights,
        }
    )
    loadings_df = pd.DataFrame({"indicator": SOIL_INDICATORS})
    for idx in range(n_components):
        loadings_df[f"PC{idx + 1}_loading"] = retained_loadings[:, idx]
        loadings_df[f"PC{idx + 1}_explained_variance"] = retained_variance[idx]
        loadings_df[f"PC{idx + 1}_eigenvalue"] = eigenvalues[idx]
    return sqi, weights_df, loadings_df


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


@dataclass
class TreeNode:
    value: float
    feature: int | None = None
    threshold: float | None = None
    left: "TreeNode | None" = None
    right: "TreeNode | None" = None


class NumpyRegressionTree:
    """A small regression-tree implementation used for reproducibility."""

    def __init__(
        self,
        *,
        max_depth: int | None = None,
        min_samples_split: int = 2,
        min_samples_leaf: int = 1,
        max_features: int | None = None,
        random_state: int | None = None,
    ) -> None:
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.max_features = max_features
        self.rng = np.random.default_rng(random_state)
        self.root: TreeNode | None = None
        self.feature_importances_: np.ndarray | None = None

    def fit(self, x: np.ndarray, y: np.ndarray) -> "NumpyRegressionTree":
        self.feature_importances_ = np.zeros(x.shape[1], dtype=float)
        self.root = self._build(x, y, depth=0)
        total = self.feature_importances_.sum()
        if total > 0:
            self.feature_importances_ /= total
        return self

    def _build(self, x: np.ndarray, y: np.ndarray, depth: int) -> TreeNode:
        value = float(np.mean(y))
        if (
            len(y) < self.min_samples_split
            or np.var(y) <= 1e-15
            or (self.max_depth is not None and depth >= self.max_depth)
        ):
            return TreeNode(value=value)

        split = self._best_split(x, y)
        if split is None:
            return TreeNode(value=value)

        feature, threshold, improvement = split
        left_mask = x[:, feature] <= threshold
        right_mask = ~left_mask
        if self.feature_importances_ is not None:
            self.feature_importances_[feature] += improvement

        return TreeNode(
            value=value,
            feature=feature,
            threshold=threshold,
            left=self._build(x[left_mask], y[left_mask], depth + 1),
            right=self._build(x[right_mask], y[right_mask], depth + 1),
        )

    def _best_split(self, x: np.ndarray, y: np.ndarray) -> tuple[int, float, float] | None:
        n_samples, n_features = x.shape
        max_features = self.max_features or n_features
        features = self.rng.choice(n_features, size=max_features, replace=False)
        parent_sse = np.sum((y - y.mean()) ** 2)
        best_feature = None
        best_threshold = None
        best_improvement = 0.0

        for feature in features:
            order = np.argsort(x[:, feature], kind="mergesort")
            x_sorted = x[order, feature]
            y_sorted = y[order]
            unique_changes = np.flatnonzero(np.diff(x_sorted) != 0)
            if len(unique_changes) == 0:
                continue

            csum = np.cumsum(y_sorted)
            csum2 = np.cumsum(y_sorted**2)
            total_sum = csum[-1]
            total_sum2 = csum2[-1]

            for idx in unique_changes:
                left_n = idx + 1
                right_n = n_samples - left_n
                if left_n < self.min_samples_leaf or right_n < self.min_samples_leaf:
                    continue
                left_sum = csum[idx]
                left_sum2 = csum2[idx]
                right_sum = total_sum - left_sum
                right_sum2 = total_sum2 - left_sum2
                left_sse = left_sum2 - (left_sum**2 / left_n)
                right_sse = right_sum2 - (right_sum**2 / right_n)
                improvement = parent_sse - left_sse - right_sse
                if improvement > best_improvement:
                    best_feature = feature
                    best_threshold = (x_sorted[idx] + x_sorted[idx + 1]) / 2.0
                    best_improvement = improvement

        if best_feature is None:
            return None
        return int(best_feature), float(best_threshold), float(best_improvement)

    def predict(self, x: np.ndarray) -> np.ndarray:
        if self.root is None:
            raise RuntimeError("Tree has not been fitted.")
        return np.array([self._predict_row(row, self.root) for row in x], dtype=float)

    def _predict_row(self, row: np.ndarray, node: TreeNode) -> float:
        while node.feature is not None and node.threshold is not None:
            if row[node.feature] <= node.threshold:
                node = node.left  # type: ignore[assignment]
            else:
                node = node.right  # type: ignore[assignment]
        return node.value


class NumpyRandomForestRegressor:
    """Small Random Forest regressor matching the manuscript calculation."""

    def __init__(
        self,
        *,
        n_estimators: int = 100,
        random_state: int = 42,
        max_features: int | None = None,
    ) -> None:
        self.n_estimators = n_estimators
        self.random_state = random_state
        self.max_features = max_features
        self.trees: list[NumpyRegressionTree] = []
        self.feature_importances_: np.ndarray | None = None

    def fit(self, x: np.ndarray, y: np.ndarray) -> "NumpyRandomForestRegressor":
        rng = np.random.default_rng(self.random_state)
        n_samples = x.shape[0]
        importances = np.zeros(x.shape[1], dtype=float)
        self.trees = []

        for _ in range(self.n_estimators):
            sample_idx = rng.integers(0, n_samples, size=n_samples)
            tree_seed = int(rng.integers(0, 2**32 - 1))
            tree = NumpyRegressionTree(
                max_depth=None,
                min_samples_split=2,
                min_samples_leaf=1,
                max_features=self.max_features or x.shape[1],
                random_state=tree_seed,
            )
            tree.fit(x[sample_idx], y[sample_idx])
            self.trees.append(tree)
            if tree.feature_importances_ is not None:
                importances += tree.feature_importances_

        total = importances.sum()
        self.feature_importances_ = importances / total if total > 0 else importances
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        if not self.trees:
            raise RuntimeError("Forest has not been fitted.")
        predictions = np.vstack([tree.predict(x) for tree in self.trees])
        return predictions.mean(axis=0)


def write_importance_svg(path: Path, importances_df: pd.DataFrame, r2: float) -> None:
    data = importances_df.sort_values("importance", ascending=True).reset_index(drop=True)
    width = 920
    row_h = 38
    margin_l = 180
    margin_r = 130
    margin_t = 82
    margin_b = 46
    height = margin_t + margin_b + row_h * len(data)
    bar_w = width - margin_l - margin_r
    max_val = float(data["importance"].max()) if len(data) else 1.0
    max_val = max(max_val, 1e-12)
    colors = ["#376795", "#4b8c61", "#d9a441", "#c95f4f"]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="28" y="34" font-family="Arial" font-size="22" font-weight="700" fill="#111">Random Forest feature importance for SQI</text>',
        f'<text x="28" y="60" font-family="Arial" font-size="13" fill="#555">n_estimators = 100; random_state = 42; training R2 = {r2:.3f}</text>',
    ]

    for idx, row in data.iterrows():
        y = margin_t + idx * row_h
        feature = escape(str(row["feature"]))
        importance = float(row["importance"])
        length = importance / max_val * bar_w
        color = colors[idx % len(colors)]
        parts.append(
            f'<text x="{margin_l - 14}" y="{y + 23}" font-family="Arial" font-size="14" text-anchor="end" fill="#222">{feature}</text>'
        )
        parts.append(
            f'<rect x="{margin_l}" y="{y + 7}" width="{length:.2f}" height="20" rx="2" fill="{color}"/>'
        )
        parts.append(
            f'<text x="{margin_l + length + 8:.2f}" y="{y + 23}" font-family="Arial" font-size="13" fill="#333">{importance:.4f}</text>'
        )

    axis_y = margin_t + row_h * len(data) + 6
    parts.append(
        f'<line x1="{margin_l}" y1="{axis_y}" x2="{margin_l + bar_w}" y2="{axis_y}" stroke="#222" stroke-width="1"/>'
    )
    parts.append(
        f'<text x="{margin_l + bar_w / 2}" y="{height - 16}" font-family="Arial" font-size="13" text-anchor="middle" fill="#333">Mean decrease in impurity (normalized)</text>'
    )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_metadata(path: Path, rows: list[tuple[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["item", "value"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydrochar-table", default="data/tableS7_hydrochar_properties.csv")
    parser.add_argument("--soil-table", default="data/tableS8_soil_quality_indicators.csv")
    parser.add_argument("--output-dir", default="outputs")
    parser.add_argument("--n-pcs", type=int, default=3)
    parser.add_argument("--n-estimators", type=int, default=100)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    hydrochar_path = Path(args.hydrochar_table)
    soil_path = Path(args.soil_table)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    hydrochar_df = pd.read_csv(hydrochar_path)
    soil_df = pd.read_csv(soil_path)
    n_records = min(len(hydrochar_df), len(soil_df))
    hydrochar_df = hydrochar_df.iloc[:n_records].reset_index(drop=True)
    soil_df = soil_df.iloc[:n_records].reset_index(drop=True)

    missing = [col for col in RF_FEATURES if col not in hydrochar_df.columns]
    if missing:
        raise ValueError(f"Hydrochar table is missing required columns: {missing}")

    y, weights_df, loadings_df = calculate_pca_sqi(soil_df, n_components=args.n_pcs)

    x_raw = hydrochar_df[RF_FEATURES].map(parse_number).to_numpy(dtype=float)

    x, impute_means = mean_impute(x_raw)
    model = NumpyRandomForestRegressor(
        n_estimators=args.n_estimators,
        random_state=args.random_state,
    )
    model.fit(x, y)
    if model.feature_importances_ is None:
        raise RuntimeError("Random Forest did not produce feature importances.")

    predictions = model.predict(x)
    train_r2 = r2_score(y, predictions)

    importance_df = pd.DataFrame(
        {"feature": RF_FEATURES, "importance": model.feature_importances_}
    ).sort_values("importance", ascending=False)

    sqi_dataset = hydrochar_df[["Reference"] + RF_FEATURES].copy()
    for feature in RF_FEATURES:
        sqi_dataset[feature] = sqi_dataset[feature].map(parse_number)
    sqi_dataset["SQI"] = y
    if "Reference" in soil_df.columns:
        sqi_dataset["Soil_reference"] = soil_df["Reference"].values
    sqi_dataset.to_csv(out_dir / "sqi_dataset_for_random_forest.csv", index=False)
    weights_df.to_csv(out_dir / "sqi_pca_weights.csv", index=False)
    loadings_df.to_csv(out_dir / "sqi_pca_loadings.csv", index=False)
    importance_df.to_csv(out_dir / "fig5c_feature_importance.csv", index=False)
    write_importance_svg(out_dir / "fig5c_feature_importance.svg", importance_df, train_r2)

    write_metadata(
        out_dir / "run_metadata.csv",
        [
            ("hydrochar_table", str(hydrochar_path)),
            ("soil_table", str(soil_path)),
            ("n_records", n_records),
            ("soil_indicators", "; ".join(SOIL_INDICATORS)),
            ("features", "; ".join(RF_FEATURES)),
            ("target", TARGET),
            ("excluded_predictor", "AR (%)"),
            ("n_pcs_for_sqi", args.n_pcs),
            ("n_estimators", args.n_estimators),
            ("random_state", args.random_state),
            ("training_r2", f"{train_r2:.6f}"),
            ("imputation", "Mean imputation for missing predictor values"),
            ("imputation_means", "; ".join(f"{f}={m:.6g}" for f, m in zip(RF_FEATURES, impute_means))),
        ],
    )

    print(f"Hydrochar table: {hydrochar_path}")
    print(f"Soil table: {soil_path}")
    print(f"Paired records: {n_records}")
    print(f"Training R2: {train_r2:.3f}")
    print("\nFeature importance:")
    for _, row in importance_df.iterrows():
        print(f"  {row['feature']}: {row['importance']:.4f}")
    print(f"\nOutputs written to: {out_dir}")


if __name__ == "__main__":
    main()

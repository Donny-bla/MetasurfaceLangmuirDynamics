import argparse
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy import stats

@dataclass
class SplitIndices:
    train: np.ndarray
    val: np.ndarray
    test: np.ndarray

def grouped_concentration_split(metadata: pd.DataFrame, args: argparse.Namespace) -> SplitIndices:
    groups = metadata["concentrationIndex"].drop_duplicates().to_numpy(copy=True)
    rng = np.random.default_rng(args.split_seed)
    rng.shuffle(groups)
    n_test = max(1, int(round(len(groups) * args.test_fraction)))
    n_val = max(1, int(round(len(groups) * args.val_fraction)))
    test_groups = set(groups[:n_test])
    val_groups = set(groups[n_test:n_test + n_val])
    group_values = metadata["concentrationIndex"].to_numpy()
    test = np.flatnonzero(np.isin(group_values, list(test_groups)))
    val = np.flatnonzero(np.isin(group_values, list(val_groups)))
    train = np.flatnonzero(~np.isin(group_values, list(test_groups | val_groups)))
    return SplitIndices(train=train, val=val, test=test)

def safe_corr(y_true: np.ndarray, y_pred: np.ndarray, method: str) -> float:
    if method == "pearson":
        return float(stats.pearsonr(y_true, y_pred).statistic)
    if method == "spearman":
        return float(stats.spearmanr(y_true, y_pred).statistic)

def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    residual = y_pred - y_true
    fold_error = np.power(10.0, np.abs(residual))
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    return {
        "rmse_log10_concentration": float(np.sqrt(np.mean(residual ** 2))),
        "mae_log10_concentration": float(np.mean(np.abs(residual))),
        "r2_log10_concentration": 1.0 - ss_res / ss_tot,
        "pearson_corr": safe_corr(y_true, y_pred, "pearson"),
        "spearman_corr": safe_corr(y_true, y_pred, "spearman"),
        "median_fold_error": float(np.median(fold_error)),
        "p90_fold_error": float(np.percentile(fold_error, 90)),
    }

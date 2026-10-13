"""Shared Monte Carlo binding simulation and physical settings loader."""
from pathlib import Path
import numpy as np
import pandas as pd

def load_dataset_config(summary_path: Path) -> dict[str, float]:
    row = pd.read_csv(summary_path).iloc[0].to_dict()
    return {
        "concentrationMinM": float(row["concentrationMinM"]),
        "concentrationMaxM": float(row["concentrationMaxM"]),
        "dtSeconds": float(row["dtSeconds"]),
        "kon_Minv_sinv": float(row["kon_Minv_sinv"]),
        "koff_sinv": float(row["koff_sinv"]),
        "nBindingSites": int(row["nBindingSites"]),
        "proteinDiameterNm": float(row["proteinDiameterNm"]),
    }

def simulate_thickness_traces(
    concentration: float,
    n_realizations: int,
    time_stop_s: int,
    output_time_s: np.ndarray,
    dt: float,
    kon: float,
    koff: float,
    n_sites: int,
    protein_diameter_nm: float,
    rng: np.random.Generator,
) -> np.ndarray:
    p_on = min(kon * concentration * dt, 1.0)
    p_off = min(koff * dt, 1.0)
    occupied = np.zeros(n_realizations, dtype=np.int32)
    output_lookup = {int(t): i for i, t in enumerate(output_time_s)}
    thickness = np.zeros((n_realizations, len(output_time_s)), dtype=np.float32)

    n_steps = int(round(time_stop_s / dt))
    for step in range(1, n_steps + 1):
        current_time = int(round(step * dt))
        n_empty = n_sites - occupied
        n_bind = rng.binomial(n_empty, p_on).astype(np.int32)
        n_unbind = rng.binomial(occupied, p_off).astype(np.int32)
        occupied += n_bind - n_unbind
        if current_time in output_lookup:
            out_idx = output_lookup[current_time]
            thickness[:, out_idx] = (np.pi / 6.0) * (occupied / n_sites) * protein_diameter_nm
    return thickness

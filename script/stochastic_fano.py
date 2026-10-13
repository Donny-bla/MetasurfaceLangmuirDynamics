from pathlib import Path

import numpy as np
import pandas as pd


NOISY_PARAMETERS = ("fanoCentreNm", "fanoGammaNm", "fanoQ", "fanoAmplitude")
PARAMETERS = (*NOISY_PARAMETERS, "fanoOffset")


def load_stochastic_fano(fano_dir: Path, wavelength_stride: int) -> dict:
    clean = pd.read_csv(fano_dir / "fano_forward_model_clean_parameters.csv")
    clean = clean.sort_values("thicknessNm")
    grid = clean["thicknessNm"].to_numpy(dtype=float)
    values = clean[list(PARAMETERS)].to_numpy(dtype=float)
    config = pd.read_csv(fano_dir / "fano_forward_model_parameter_noise_config.csv").set_index("parameterName")
    std = config.loc[list(NOISY_PARAMETERS), "noiseStd"].to_numpy(dtype=float)
    wavelengths = pd.read_csv(
        fano_dir / "fano_forward_model_clean_spectrum_matrix.csv", usecols=["wavelengthNm"]
    )["wavelengthNm"].to_numpy(dtype=float)[::wavelength_stride]
    coefficients = np.polynomial.polynomial.polyfit(grid, values, 2)
    return dict(wavelength_nm=wavelengths, thickness_nm=grid, parameters=values,
                coefficients=coefficients, parameter_noise_std=std)


def sample_parameters(thickness: np.ndarray, fano: dict, rng: np.random.Generator) -> np.ndarray:
    thickness = np.asarray(thickness, dtype=float)
    grid = fano["thickness_nm"]
    clipped = np.clip(thickness, grid[0], grid[-1])
    params = np.moveaxis(np.polynomial.polynomial.polyval(clipped, fano["coefficients"]), 0, -1)
    params[..., :4] += rng.normal(size=(*thickness.shape, 4)) * fano["parameter_noise_std"]
    params[..., 1] = np.maximum(params[..., 1], np.finfo(float).eps)
    params[..., 3] = np.maximum(params[..., 3], 0.0)
    return params


def sample_spectral_images(thickness: np.ndarray, fano: dict, rng: np.random.Generator) -> np.ndarray:
    params = sample_parameters(thickness, fano, rng)
    centre, gamma, q, amplitude, offset = [params[..., i, None] for i in range(5)]
    epsilon = (fano["wavelength_nm"] - centre) / gamma
    spectra = offset - amplitude * (q + epsilon) ** 2 / ((1 + q ** 2) * (1 + epsilon ** 2))
    return spectra.transpose(0, 2, 1).astype(np.float32)


def noise_metadata(fano: dict) -> dict:
    return {"distribution": "independent Gaussian before gamma/amplitude constraints",
            "sampling": "fresh per record and time point; shared across wavelengths",
            "noise_std": dict(zip(NOISY_PARAMETERS, fano["parameter_noise_std"].tolist())),
            "offset_noise_std": 0.0,
            "clean_parameter_evaluation": "quadratic"}

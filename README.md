# MetasurfaceLangmuirDynamics

Supporting repository for the manuscript on time-resolved full-spectrum decoding for quasi-BIC metasurface biosensing.

This project models adsorption-driven spectral changes in a quasi-bound-state-in-the-continuum (quasi-BIC) metasurface and evaluates whether full wavelength-time spectra contain more useful information for concentration estimation than conventional resonance tracking alone.

## Overview

The workflow combines:

- metasurface optical calibration using Fano-line fitting,
- stochastic Langmuir adsorption dynamics for protein binding,
- noisy time-resolved transmission spectra generated from the binding trajectory,
- concentration regression using peak-based and full-spectrum decoders,
- comparison against a structured Monte Carlo Fisher-information reference.

The repository is intended as a reproducibility and supplementary resource for the manuscript rather than as a standalone software package.

## Project goal

The core question is whether the entire wavelength-time response of a biosensing metasurface can provide more accurate and robust protein concentration estimation than resonance-peak tracking alone.

The supporting analysis includes:

- comparison of metasurface geometries for bulk vs. surface sensitivity,
- calibration of the Fano response versus effective protein-layer thickness,
- generation of a synthetic AFP concentration dataset under stochastic Langmuir kinetics,
- benchmarking of four decoders: peak+Langmuir, peak+SVR, bilinear regression, and a 2D CNN,
- concentration-resolved CRB analysis.

## Repository structure

```text
MetasurfaceLangmuirDynamics/
├── README.md
├── inputs/
│   ├── fit_coated/
│   ├── fit_zero/
│   └── raw/
├── script/
│   ├── __init__.py
│   ├── analyzeFanoThicknessSpectra.m
│   ├── binding_helpers.py
│   ├── bilinear_spectral_regressor.ipynb
│   ├── build_double_ellipse_fano_forward_model.ipynb
│   ├── build_langmuir_spectral_images.ipynb
│   ├── concentration_split.py
│   ├── estimate_langmuir_crb_structured_mc.ipynb
│   ├── fitITFFano.m
│   ├── fit_current_fano.m
│   ├── generate_fixed_rate_langmuir_dataset.ipynb
│   ├── peak_trace_svr_regressor.ipynb
│   ├── peak_tracking_langmuir_fit.ipynb
│   ├── simulate_current_langmuir.m
│   ├── stochastic_fano.py
│   └── train_langmuir_2d_cnn_regressor.ipynb
└── ...
```

## Typical workflow

### 1) Fit the optical response model

The optics and spectral fitting pipeline is implemented primarily in MATLAB scripts and notebooks:

- `script/fitITFFano.m`
- `script/analyzeFanoThicknessSpectra.m`
- `script/fit_current_fano.m`
- `script/build_double_ellipse_fano_forward_model.ipynb`

These steps fit the quasi-BIC transmission spectra to a Fano model and build the dependence of the fitted parameters on the effective protein-layer thickness.

### 2) Generate Langmuir-based spectral trajectories

The binding and data-generation flow is centered on:

- `script/simulate_current_langmuir.m`
- `script/stochastic_fano.py`
- `script/generate_fixed_rate_langmuir_dataset.ipynb`
- `script/build_langmuir_spectral_images.ipynb`

This produces noisy wavelength-time spectral images for a range of concentrations, following stochastic adsorption/desorption dynamics and the Fano forward model.

### 3) Train and compare concentration decoders

The machine-learning experiments are in:

- `script/peak_tracking_langmuir_fit.ipynb`
- `script/peak_trace_svr_regressor.ipynb`
- `script/bilinear_spectral_regressor.ipynb`
- `script/train_langmuir_2d_cnn_regressor.ipynb`

These notebooks benchmark the performance of:

- peak tracking + Langmuir fitting,
- peak tracking + RBF-SVR,
- rank-2 bilinear regression on full spectral images,
- 2D CNN on wavelength-time data.

### 4) Evaluate information content and estimator bounds

The structured Monte Carlo CRB analysis is in:

- `script/estimate_langmuir_crb_structured_mc.ipynb`
- `script/concentration_split.py`

This compares empirical estimator RMSE with a local Fisher-information-based lower bound for log-concentration estimation.

## Data organization

The repository includes calibration and experimental input folders under `inputs/`:

- `inputs/raw/`: raw spectral and/or measurement input data
- `inputs/fit_zero/`: zero-layer or baseline fitting data
- `inputs/fit_coated/`: fitted or coated-state spectral data used for calibration

Exact file names and formats depend on the specific experiment or computational run and may need to be aligned with the scripts before execution.

## Dependencies

The project uses a mix of MATLAB and Python workflows.

Python dependencies commonly required include:

- Python 3.9+
- NumPy
- SciPy
- Matplotlib
- pandas
- scikit-learn
- PyTorch

MATLAB dependencies include standard signal-processing and optimization tools used in the Fano-fitting scripts.

## Recommended usage

1. Open the Fano-fitting scripts to calibrate the transmission model against thickness-dependent spectra.
2. Use the Langmuir and stochastic Fano scripts to generate concentration-dependent time-resolved spectral datasets.
3. Train the decoders using the matching notebooks.
4. Run the CRB notebook to compare model error against the information-theoretic lower bound.
5. Use the generated results and figures to support the manuscript analysis.

## Notes for reproducibility

- This repository is a supporting codebase for a simulation study and should be treated as a research supplement.
- The scripts are organized as an exploratory analysis workflow and may require minor path or environment adjustments depending on local machine setup.
- The concentration range, binding kinetics, and noise assumptions used in the study are defined in the manuscript and in the notebook code.

## Citation

If you use this repository in a project or publication, please cite the associated manuscript and reference the repository URL.

## Contact

For questions or collaboration, please contact the corresponding author or repository maintainer listed in the manuscript.


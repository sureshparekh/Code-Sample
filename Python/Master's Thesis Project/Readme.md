# Extended X-ray Jets in Radio-loud Quasars: A Morphological Study

`Bayesian Analysis`, `Markov Chain Monte Carlo Simulation`, `LIRA`, `pylira`

The energy produced by the supermassive black hole (SMBH) at the centre of active galactic nuclei is transported over vast distances (>100 kpcs) by jets. In the early universe, the effects of jets on the surrounding environment play a role in the development and evolution of structures. When seen at modest angles to the line of sight, radio-loud quasars’ radiation from the innermost jets (parsec scales or smaller) can be Doppler boosted.

Relativistic Jets can therefore give useful observational information of the state of SMBH activity. Over their lifetime, large-scale X-ray jets record the history of SMBH activity at distances of up to hundreds of kiloparsecs. The Chandra X-ray Observatory, due to its excellent imaging capabilities, has successfully resolved jets at X-ray energy band (0.5 - 7 keV).

Even while the number of these X-ray jets has significantly increased since Chandra’s launch in 1999, it is still relatively minor compared to the total number of recognised quasars. Only a small number of the about 100 large-size X-ray jets discovered so far have high-quality X-ray morphology data. Detecting X-ray jets is intrinsically difficult because of the low number of X-ray photon counts obtained from jets compared to their respective quasar cores. Therefore identifying and characterizing extended X-ray jets is of at most importance to understand their morphology, particle acceleration process and radiation mechanisms responsible for the observation of these enigmatic structures.

To achieve this objective, we have used the publicly available tool, multi- scale Bayesian method known as Low count Image Reconstruction and Analysis (LIRA). The algorithm models the residual as a multi-scale component and generates a series of images that capture the emission that may be present more than the baseline model. We can then compute a p-value by generating a series of **Markov Chain Monte Carlo simulations** of images under the baseline model and fitting each of these simulated images using LIRA.

## What is LIRA?

The Low-counts Image Reconstruction and Analysis (LIRA) is a robust statistical tool that is designed to estimate the statistical significance of emission from low-count images (e.g., shallow *Chandra* observations) while accounting for emission from surrounding bright objects and Poisson background. See its github [page](https://github.com/astrostat/LIRA) for more details and references.

## `pylira`

`pylira` is a Python package for Bayesian low-counts image reconstruction and analysis. See its github [page](https://github.com/astrostat/pylira) for more details and references.

### Installation

---

Pylira requires a working **R** installation to build against. On Mac-OS you can typically use:

```
brew install r
```

On Linux based systems the following should work:

```
sudo apt-get install r-base-dev r-base r-mathlib
```

Once `R` is installed you can install Pylira using:

```
pip install pylira
```

## Prerequisites

* Latest version of CIAO (conda installation [instructions](https://cxc.harvard.edu/ciao/download/conda.html))
* Latest version of Python compatable with CIAO
* Latest version of SAOTrace, `sherpa`, `marx`, and `astropy`


## Input preparation

Python helper scripts to prepare inputs are located in the folder. These helper scripts are tailored for use with *Chandra* observations.

 After configuring ``lira_input_prep.yaml`` in the folder with inputs (event files, region files), exceute ``create_inputs_for_lira.py``. This will:

* Create a binned image from the events file
* Extract the spectrum of the core
* Fit an absorbed powerlaw model to the spectrum
* Create a PSF of the core with matched pixel boundaries
* Create a baseline image and the specified number of replicates.

## Execution

The following input files are required before running LIRA:

* Input observation (size should be a power of two)
* PSF image
* Baseline image (usually point source+flat background; same size as the input image)
* Baseline replicates (*faked* images of the baseline model; same size as the input image)
* ROI/region files (CIAO format) to test for significance
* ``config.yaml``
* Setup the `run_lira_final.py` as per the inputs and run.
* Use the `plot_lira.py` to plot the data files generated by `run_lira_final.py`.

To understand the importance and usage of this pipeline, refer to the `presentation.pdf`. You can find the detailed report in `Masters_Thesis.pdf`.

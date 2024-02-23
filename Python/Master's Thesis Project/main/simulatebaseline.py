from prepare_spectra import prepare_spectra
import yaml
from rebin_img import rebin_img
from rebin_psf import rebin_psf
from ciao_contrib.runtool import *
from sherpa.astro.ui import *
from sherpa_contrib.chart import *

# from sim_image import sim_image
from saotrace_helpers import (
    run_sao_raytrace,
    run_sao_raytrace_parallel,
    run_sao_raytrace_parallel_flux,
)
import logging
import numpy as np
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import shutil


# slightly modified version of kmc's code from astrostat/LIRA

def simulate_null_images(infile, psffile, num_sims, no_core, mcmciter=5000):
    print("Creating the null file")
    clean()
    set_stat("cstat")
    set_method("simplex")
    load_image(infile)
    load_psf("mypsf", psffile)
    set_psf(mypsf)

    if no_core:
        set_model(const2d.c0)
        set_par(c0.c0, min=0)
    else:
        set_model(gauss2d.q1 + const2d.c0)
        set_par(c0.c0, min=0)
        # set_par(q1.fwhm,max=0.5)
    guess(q1)
    fit()
    results = get_fit_results()
    save("core_source_fit.save", clobber=True)
    save_source("null_q1_c1.fits", clobber=True)
    covar()

    if no_core:
        for i in range(num_sims):
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        clean()
        return 0

    normgauss1d.g1
    g1.pos = q1.fwhm
    g1.fwhm = get_covar_results().parmaxes[0]

    # check if there is a valid upper bound.
    print(get_covar_results())
    if (
        get_covar_results().parmaxes[0] is None
        or get_covar_results().parmins[1] is None
        or get_covar_results().parmins[0] is None
    ):
        for i in range(num_sims):
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        clean()
        return 0
    # if not go for the regular
    set_prior(q1.fwhm, g1)
    set_sampler_opt("defaultprior", False)
    set_sampler_opt("priorshape", [True, False, False, False, False])
    set_sampler_opt("originalscale", [True, True, True, True, True])
    if mcmciter < num_sims * 100:
        mcmciter = num_sims * 100

    # the following code throws an error sometimes #bug
    try:
        stats, accept, params = get_draws(1, niter=mcmciter)
    except:
        params = [np.repeat(q1.fwhm.val, mcmciter)]

    # print('Simulating the null files')
    for i in range(num_sims):
        set_par(q1.fwhm, params[0][(i + 1) * 100 - 1])
        fake()
        save_image("sim_null_{}.fits".format(i), clobber=True)
    save_all(outfile="lira_input_baseline_sim.log", clobber=True)
    clean()

simulate_null_images(
            "img_64x64_0.492.fits",
            "psf_demo.fits",
            50,
            True,
        )
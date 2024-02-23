# import os
# import tempfile
from copy import deepcopy
# from pathlib import Path
import numpy as np
# from scipy.ndimage import labeled_comprehension
# from astropy.table import Table
# from . import image_analysis

#inputs
# with open('/home/suresh/Chandra/PKS0605018/Codes/result_obs_cen_itr1_region_0.txt', 'r') as file:
xi_obs = np.loadtxt('/home/suresh/Chandra/PKS0605018/Codes/result_obs_jet_itr1_region_2.txt')
print(xi_obs)

# with open('/home/suresh/Chandra/PKS0605018/Codes/mean_pvalues_sim_cen_region_0.txt', 'r') as file:
xi_repl = np.loadtxt('/home/suresh/Chandra/PKS0605018/Codes/mean_pvalues_sim_jet_region_2.txt')
print(xi_repl)


def _plot_xi(xi_dist, ax, ls="--", c="gray", tol=1e-10, label=None):
    from scipy import stats

    xi_dist_c = deepcopy(xi_dist)
    xi_dist_c[xi_dist_c <= tol] = tol

    xi_dist_c = np.log10(xi_dist_c)

    kernel = stats.gaussian_kde(xi_dist_c)
    eval_points = np.linspace(np.min(xi_dist_c), 0, 100)
    kde = kernel(eval_points)

    ax.plot(eval_points, kde, ls=ls, c=c, label=label)

def plot_xi_dist(xi_obs, xi_repl, figsize=(8, 5)):
    """
    Plot the posterior distributions of xi for a region

    Parameters
    ----------
    xi_obs : `~numpy.ndarray`
        Posterior distribution of xi for the observation
    xi_repl : `~numpy.ndarray`
        Posterior distribution of xi for all the replicates
    region_id : int
        Integer representing the region
    figsize : tuple
        Figure size
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    n_replicates = int(xi_repl.shape[0] / xi_obs.shape[0])
    n_iters = xi_obs.shape[0]

    # plot the replicate distribution
    for i in range(0, n_replicates * n_iters, n_iters):
        _plot_xi(xi_repl[i : i + n_iters], ax)

    # plot the mean distribution
    _plot_xi(
        xi_repl, ax, ls="-", c="black", label="Mean null distribution"
    )

    # plot the observed distribution
    _plot_xi(
        xi_obs, ax, ls="-", c="blue", label="Best fit distribution"
    )

    ax.set_xlabel(r"Posterior distribution (log$_{10}\xi$)")
    ax.set_ylabel("Density")

    ax.set_title(f"Region: 0.0")
    plt.legend()
    plt.savefig('demo.png')
    return ax

plot_xi_dist(
    xi_obs=xi_obs,
    xi_repl=xi_repl,
    # region_id="1.0",
    # plotname = 'for_central_region.png'
);

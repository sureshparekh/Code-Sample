from astropy.io import fits
import copy
import numpy as np
import matplotlib.pyplot as plt
from regions import CirclePixelRegion
from astropy.nddata.utils import Cutout2D
from astropy.io import fits
from multiprocessing import Pool
from copy import deepcopy
from pylira.data import point_source_gauss_psf, gauss_and_point_sources_gauss_psf
from pylira.utils.plot import plot_example_dataset, plot_hyperpriors_lira
from pylira import LIRADeconvolver, LIRADeconvolverResult, LIRASignificanceEstimator
from astropy.nddata.utils import Cutout2D
from regions import CirclePixelRegion, PixCoord
# from astropy.io import fits
random_state = np.random.RandomState(148)



psf = fits.getdata("/home/suresh/Chandra/PKS0605018/12056_repro/core_psf.fits")
background = fits.getdata("/home/suresh/Chandra/PKS0605018/12056_repro/background.fits")
exp_64 = fits.getdata("/home/suresh/Chandra/PKS0605018/12056_repro/expmap/expmap_64.expmap")

xi_sim_cen = []
xi_sim_cen_0 = []
xi_sim_cen_1 = []
xi_sim_cen_2 = []
xi_sim_jet = []
xi_sim_jet_0 = []
xi_sim_jet_1 = []
xi_sim_jet_2 = []


#first for the central region for simulated images
for s in range (50):
    print("running for simulated central region iteration:{0}".format(s))
    inpfile = "/home/suresh/Chandra/PKS0605018/12056_repro/sim_null_{0}.fits".format(s)
    counts = fits.getdata(inpfile)

    #input data file
    data = {
    "counts": counts.astype(np.float32),
    "psf": psf.astype(np.float32),
    "background": background.astype(np.float32),
    "exposure": exp_64.astype(np.float32),
    }

    alpha_init = 0.1 * np.ones(np.log2(data["counts"].shape[0]).astype(int))

    dec = LIRADeconvolver(
        alpha_init=alpha_init,
        n_iter_max=5000,
        n_burn_in=4000,
        fit_background_scale=False,
        random_state=random_state
    )

    data["flux_init"] = random_state.gamma(30, size=data["counts"].shape)

    print('running first simulation')
    result = dec.run(data)
    
    #Significance estimation / null hyphothesis simulation
    # number of simulations
    n_iter = 10

    # reduce number of iterations for LIRA
    dec.n_iter_max =1000
    datasets = []
    data_null = copy.deepcopy(data)

    for idx in range(n_iter):
        data_null["counts"] = random_state.poisson(data["background"])
        data_null["flux_init"] = random_state.gamma(10, size=data["counts"].shape)
        datasets.append(data_null)
    
    #time
    with Pool(processes=8) as pool:
        results= pool.map(dec.run, datasets)

    # Region of interest for the source 
    region_src_cen = CirclePixelRegion(
    center=PixCoord(33, 33),
    radius=5
    )

    # Some control region 
    region_bkg = CirclePixelRegion(
    center=PixCoord(39, 55),
    radius=5
    )
    #defining region of interest
    labels = np.zeros((64, 64))
   
    for idx, region in enumerate([region_src_cen, region_bkg]):
        mask = region.to_mask()
        labels += (idx + 1) * mask.to_image(shape=labels.shape)
    
    #estimation and visualization of Xi. First we instantiate the LIRASignificanceEstimator with the result obtained before, the results from the null hypothesis simulations and the label image
    est = LIRASignificanceEstimator(
    result_observed_im=result,
    result_replicates=results,
    labels_im=labels
    )

    #Now we first estimate the p-values for the regions defined by the labels.
    result_p_values_cen = est.estimate_p_values(
        data=data, gamma=0.005
    )
    
    p_dis_c_0=result_p_values_cen[2]['0.0']
    p_dis_c_1=result_p_values_cen[2]['1.0']
    p_dis_c_2=result_p_values_cen[2]['2.0']
    #for k in range(len(p_dis_c)):
    #    xi_sim+=[p_dis_c[k]]
    # np.savetxt("result_centralsim_itr_{0}.txt".format(s), result_p_values_cen[2])    
    np.savetxt("result_sim_cen_itr1_region_0_{0}.txt".format(s), p_dis_c_0) 
    np.savetxt("result_sim_cen_itr1_region_1_{0}.txt".format(s), p_dis_c_1) 
    np.savetxt("result_sim_cen_itr1_region_2_{0}.txt".format(s), p_dis_c_2) 
    xi_sim_cen_0.append(p_dis_c_0)
    xi_sim_cen_1.append(p_dis_c_1)
    xi_sim_cen_2.append(p_dis_c_2)
    # del xi_sim

    # xi_sim_cen_mean = np.mean(xi_sim_cen, axis=0)
    # np.savetxt("xi_sim_cen_mean.txt", xi_sim_cen_mean)

    # xi_sim_cen.append(result_p_values_cen[2])
    # np.savetxt("result_centralsim_itr1_{0}.txt".format(s), xi_sim_cen) 
        # Region of interest for the source 
    region_src_jet = CirclePixelRegion(
    center=PixCoord(16, 31),
    radius=5
    )

    # Some control region 
    region_bkg = CirclePixelRegion(
    center=PixCoord(39, 55),
    radius=5
    )
    #defining region of interest
    labels = np.zeros((64, 64))

    for idx, region in enumerate([region_src_jet, region_bkg]):
        mask = region.to_mask()
        labels += (idx + 1) * mask.to_image(shape=labels.shape)

    #estimation and visualization of Xi. First we instantiate the LIRASignificanceEstimator with the result obtained before, the results from the null hypothesis simulations and the label image
    est_jet = LIRASignificanceEstimator(
    result_observed_im=result,
    result_replicates=results,
    labels_im=labels
    )

    #Now we first estimate the p-values for the regions defined by the labels.
    result_p_values_jet = est_jet.estimate_p_values(
        data=data, gamma=0.005
    )
    # for k in range(len(result_p_values_jet[2])):
        # xi_sim+=[result_p_values_jet[2][k]]
    # np.savetxt("result_centralsim_itr_{0}.txt".format(s), result_p_values_jet[2])    
    # np.savetxt("result_centralsim_itr1_{0}.txt".format(s), xi_sim) 
    
    p_dis_j_0=result_p_values_jet[2]['0.0']
    p_dis_j_1=result_p_values_jet[2]['1.0']
    p_dis_j_2=result_p_values_jet[2]['2.0']
    #for k in range(len(p_dis_c)):
    #    xi_sim+=[p_dis_c[k]]
    # np.savetxt("result_centralsim_itr_{0}.txt".format(s), result_p_values_cen[2])    
    np.savetxt("result_sim_jet_itr1_region_0_{0}.txt".format(s), p_dis_j_0) 
    np.savetxt("result_sim_jet_itr1_region_1_{0}.txt".format(s), p_dis_j_1) 
    np.savetxt("result_sim_jet_itr1_region_2_{0}.txt".format(s), p_dis_j_2) 
    xi_sim_jet_0.append(p_dis_j_0)
    xi_sim_jet_1.append(p_dis_j_1)
    xi_sim_jet_2.append(p_dis_j_2)



    # xi_sim_jet.append(result_p_values_jet[2])
    # np.savetxt("result_centralsim_itr1_{0}.txt".format(s), xi_sim_jet) 

# np.savetxt("appended_pvalues_region_0.txt", xi_sim_cen_0) 
# np.savetxt("appended_pvalues_region_1.txt", xi_sim_cen_1) 
# np.savetxt("appended_pvalues_region_2.txt", xi_sim_cen_2) 
# np.savetxt("appended_pvalues_region_0.txt", xi_sim_jet_0) 
# np.savetxt("appended_pvalues_region_1.txt", xi_sim_jet_1) 
# np.savetxt("appended_pvalues_region_2.txt", xi_sim_jet_2) 
xi_sim_cen_mean_0 = np.mean(xi_sim_cen_0, axis=0)
xi_sim_cen_mean_1 = np.mean(xi_sim_cen_1, axis=0)
xi_sim_cen_mean_2 = np.mean(xi_sim_cen_2, axis=0)

xi_sim_jet_mean_0 = np.mean(xi_sim_cen_0, axis=0)
xi_sim_jet_mean_1 = np.mean(xi_sim_cen_1, axis=0)
xi_sim_jet_mean_2 = np.mean(xi_sim_cen_2, axis=0)

np.savetxt("mean_pvalues_sim_cen_region_0.txt", xi_sim_cen_mean_0) 
np.savetxt("mean_pvalues_sim_cen_region_1.txt", xi_sim_cen_mean_1) 
np.savetxt("mean_pvalues_sim_cen_region_2.txt", xi_sim_cen_mean_2) 
np.savetxt("mean_pvalues_sim_jet_region_0.txt", xi_sim_jet_mean_0) 
np.savetxt("mean_pvalues_sim_jet_region_1.txt", xi_sim_jet_mean_1) 
np.savetxt("mean_pvalues_sim_jet_region_2.txt", xi_sim_jet_mean_2) 


#for obs is below
xi_obs_cen_0=[]
xi_obs_cen_1=[]
xi_obs_cen_2=[]

xi_obs_jet_0=[]
xi_obs_jet_1=[]
xi_obs_jet_2=[]

counts = fits.getdata("/home/suresh/Chandra/PKS0605018/12056_repro/img_64x64_0.5.fits") # this should be real observation.

#input data file
data = {
"counts": counts.astype(np.float32),
"psf": psf.astype(np.float32),
"background": background.astype(np.float32),
"exposure": exp_64.astype(np.float32),
}

alpha_init = 0.1 * np.ones(np.log2(data["counts"].shape[0]).astype(int))

dec = LIRADeconvolver(
    alpha_init=alpha_init,
    n_iter_max=5000,
    n_burn_in=4000,
    fit_background_scale=False,
    random_state=random_state
)

data["flux_init"] = random_state.gamma(30, size=data["counts"].shape)

print("running for observation central region")
result = dec.run(data)

#Significance estimation / null hyphothesis simulation
# number of simulations
n_iter = 10

# reduce number of iterations for LIRA
dec.n_iter_max =1000
datasets = []
data_null = copy.deepcopy(data)

for idx in range(n_iter):
    data_null["counts"] = random_state.poisson(data["background"])
    data_null["flux_init"] = random_state.gamma(10, size=data["counts"].shape)
    datasets.append(data_null)


#time
with Pool(processes=8) as pool:
    results= pool.map(dec.run, datasets)

# Region of interest for the source 
region_src_cen = CirclePixelRegion(
center=PixCoord(33, 33),
radius=5
)

# Some control region 
region_bkg = CirclePixelRegion(
center=PixCoord(39, 55),
radius=5
)
#defining region of interest
labels = np.zeros((64, 64))

for idx, region in enumerate([region_src_cen, region_bkg]):
    mask = region.to_mask()
    labels += (idx + 1) * mask.to_image(shape=labels.shape)

#estimation and visualization of Xi. First we instantiate the LIRASignificanceEstimator with the result obtained before, the results from the null hypothesis simulations and the label image
est_cen = LIRASignificanceEstimator(
result_observed_im=result,
result_replicates=results,
labels_im=labels
)

#Now we first estimate the p-values for the regions defined by the labels.
result_p_values_cen = est.estimate_p_values(
    data=data, gamma=0.005
)
# for k in range(len(result_p_values_cen[2])):
    # xi_sim+=[result_p_values_cen[2][k]]
# np.savetxt("result_centralobs_itr.txt", result_p_values_cen[2])    
# np.savetxt("result_centralobs_itr1.txt", xi_sim) 
# p_disc_c =result_p_values_cen[2] 
p_dis_c_0=result_p_values_cen[2]['0.0']
p_dis_c_1=result_p_values_cen[2]['1.0']
p_dis_c_2=result_p_values_cen[2]['2.0']
np.savetxt("result_obs_cen_itr1_region_0.txt", p_dis_c_0) 
np.savetxt("result_obs_cen_itr1_region_1.txt", p_dis_c_1) 
np.savetxt("result_obs_cen_itr1_region_2.txt", p_dis_c_2) 

# xi_obs_cen.append(p_disc_c)
# xi_obs_cen_0.append(p_dis_c_0)
# xi_obs_cen_1.append(p_dis_c_1)
# xi_obs_cen_2.append(p_dis_c_2)
# np.savetxt("mean_obs_cen_itr1_region_0.txt", xi_obs_cen_0)
# np.savetxt("mean_obs_cen_itr1_region_1.txt", xi_obs_cen_1)
# np.savetxt("mean_obs_cen_itr1_region_2.txt", xi_obs_cen_2) 
# Region of interest for the source 
region_src_jet = CirclePixelRegion(
center=PixCoord(16, 31),
radius=5
)

# Some control region 
region_bkg = CirclePixelRegion(
center=PixCoord(39, 55),
radius=5
)
#defining region of interest
labels = np.zeros((64, 64))

for idx, region in enumerate([region_src_jet, region_bkg]):
    mask = region.to_mask()
    labels += (idx + 1) * mask.to_image(shape=labels.shape)

#estimation and visualization of Xi. First we instantiate the LIRASignificanceEstimator with the result obtained before, the results from the null hypothesis simulations and the label image
est_jet = LIRASignificanceEstimator(
result_observed_im=result,
result_replicates=results,
labels_im=labels
)

#Now we first estimate the p-values for the regions defined by the labels.
result_p_values_jet = est_jet.estimate_p_values(
    data=data, gamma=0.005
)

# p_disc_j = result_p_values_jet[2]
p_dis_j_0=result_p_values_jet[2]['0.0']
p_dis_j_1=result_p_values_jet[2]['1.0']
p_dis_j_2=result_p_values_jet[2]['2.0']

np.savetxt("result_obs_jet_itr1_region_0.txt", p_dis_j_0) 
np.savetxt("result_obs_jet_itr1_region_1.txt", p_dis_j_1) 
np.savetxt("result_obs_jet_itr1_region_2.txt", p_dis_j_2) 

# xi_obs_jet.append(p_disc_j)

# xi_obs_jet_0.append(p_dis_j_0)
# xi_obs_jet_1.append(p_dis_j_1)
# xi_obs_jet_2.append(p_dis_j_2)
# np.savetxt("mean_obs_jet_itr1_region_0.txt", xi_obs_jet_0)
# np.savetxt("mean_obs_jet_itr1_region_1.txt", xi_obs_jet_1)
# np.savetxt("mean_obs_jet_itr1_region_2.txt", xi_obs_jet_2) 


'''
#Xi plot for central region
ax = plot_xi_dist(
    xi_obs=xi_sim_cen_mean_0,
    xi_repl=p_dis_j_0,
    # region_id="0.0",
    # plotname = 'for_central_region.png'
);



ax = est_cen.plot_xi_dist(
    xi_obs=xi_obs_cen,
    xi_repl=xi_sim_cen_mean,
    region_id="1.0",
    # plotname = 'for_central_region.png'
);

ax = est_cen.plot_xi_dist(
    xi_obs=xi_obs_cen,
    xi_repl=xi_sim_cen_mean,
    region_id="2.0",
    # plotname = 'for_central_region.png'
);

#Xi plot for jet region
ax = est_jet.plot_xi_dist(
    xi_obs=xi_obs_cen,
    xi_repl=xi_sim_cen_mean,
    region_id="0.0",
    # plotname = 'for_jet_region.png'
);

#Xi plot for jet region
ax = est_jet.plot_xi_dist(
    xi_obs=xi_obs_cen,
    xi_repl=xi_sim_cen_mean,
    region_id="1.0",
    # plotname = 'for_jet_region.png'
);
#Xi plot for jet region
ax = est_jet.plot_xi_dist(
    xi_obs=xi_obs_cen,
    xi_repl=xi_sim_cen_mean,
    region_id="2.0",
    # plotname = 'for_jet_region.png'
);

'''





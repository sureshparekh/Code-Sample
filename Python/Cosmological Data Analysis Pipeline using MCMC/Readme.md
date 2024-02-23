# Observational Cosmology Project

This project focuses on utilizing observational data sets including Hz, BAO, and Pantheon to estimate parameter values of cosmological models. Subsequently, the age of the universe, the Hubble constant, and cosmography parameters are formulated based on the estimated parameters. The main methodology employed in this project is the Markov Chain Monte Carlo (MCMC) method for parameter estimation. You can refer to our publication [here](https://arxiv.org/abs/2310.18665) for more detail about the collaboration and the project. This pipeline is flexible and easy to use with any dataset and all types of cosmological models.

## Flow of the Project

1. **Data Preparation**:

   - Import necessary libraries including NumPy, Pandas, and Matplotlib.
   - Load the observational data sets (Hz, BAO, Pantheon).
2. **Model Definition**:

   - Define the cosmological model equation based on the parameters to be estimated.
3. **Likelihood Calculation**:

   - Define the log-likelihood function using the model equation to evaluate how well the model fits the data.
4. **Maximum Likelihood Estimation (MLE)**:

   - Use optimization techniques (e.g., minimize negative log-likelihood) to find the parameter values that maximize the likelihood.
5. **Prior Specification**:

   - Define the prior distributions for the parameters to incorporate prior knowledge or assumptions.
6. **Posterior Calculation**:

   - Combine the likelihood and prior to calculate the posterior distribution of the parameters.
7. **MCMC Sampling**:

   - Implement the Metropolis-Hastings algorithm to sample from the posterior distribution.
8. **Parameter Estimation**:

   - Extract parameter estimates from the MCMC chain, discarding burn-in samples and thinning for efficient sampling.
9. **Cosmological Inferences**:

   - Calculate the Hubble constant (Ho) based on the estimated parameters
   - Formulate the age of the universe, Hubble constant, and cosmography parameters using derived quantities.
10. **Cosmography**:

    - Use `plots.ipynb` file to calculate and plot other cosmographical parameters.

## About the Developer

This project is developed by Suresh Parekh. For inquiries, feedback, or collaborations, please contact [thesureshparekh@gmail.com](mailto:thesureshparekh@gmail.com).

---

*This project is hosted on a local environment and maintained by Suresh Parekh.*

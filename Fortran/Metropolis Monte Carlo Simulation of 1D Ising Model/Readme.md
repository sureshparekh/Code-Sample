# Solving 1-D Ising Model using Metropolis Monte Carlo Method



The Ising model can often be difficult to evaluate numerically if there are many states in the system. Consider an Ising model with

*L* = |Λ|: the total number of sites on the lattice,

σ~*j*~ ∈ {−1, +1}: an individual spin site on the lattice, *j* = 1, ...,  *L,*

S ∈ {−1, +1} ^*L*^ : state of the system.

Since every spin site has ±1 spin, there are *2*^*L*^ different states that are possible. This motivates the reason for the Ising model to be simulated using Monte Carlo Methods.

The Hamilton that is commonly used to represent the energy of the model when using Monte Carlo methods is

![{\displaystyle H(\sigma )=-J\sum _{\langle i~j\rangle }\sigma _{i}\sigma _{j}-h\sum _{j}\sigma _{j}.}](https://wikimedia.org/api/rest_v1/media/math/render/svg/49a000075b3808316f1976f25c235377a083841f)

Furthermore, the Hamiltonian is further simplified by assuming zero external field  *h* , since many questions that are posed to be solved using the model can be answered in absence of an external field. This leads us to the following energy equation for state σ:

![{\displaystyle H(\sigma )=-J\sum _{\langle i~j\rangle }\sigma _{i}\sigma _{j}.}](https://wikimedia.org/api/rest_v1/media/math/render/svg/e980ea997317a8166cdf80b3e27bce8eca5a57d7)

Given this Hamiltonian, quantities of interest such as the specific heat or the magnetization of the magnet at a given temperature can be calculated.

### Metropolis algorithm

The Metropolis–Hastings algorithm is the most commonly used Monte Carlo algorithm to calculate Ising model estimations. The algorithm first chooses *selection probabilities**g* (μ, ν), which represent the probability that state ν is selected by the algorithm out of all states, given that one is in state μ. It then uses acceptance probabilities  *A* (μ, ν) so that detailed balance is satisfied. If the new state ν is accepted, then we move to that state and repeat with selecting a new state and deciding to accept it. If ν is not accepted then we stay in μ. This process is repeated until some stopping criterion is met, which for the Ising model is often when the lattice becomes ferromagnetic, meaning all of the sites point in the same direction.

When implementing the algorithm, one must ensure that  *g* (μ, ν) is selected such that ergodicity is met. In thermal equilibrium a system's energy only fluctuates within a small range. This is the motivation behind the concept of  **single-spin-flip dynamics** , which states that in each transition, we will only change one of the spin sites on the lattice. Furthermore, by using single- spin-flip dynamics, one can get from any state to any other state by flipping each site that differs between the two states one at a time.

The maximum amount of change between the energy of the present state, *H*~μ~ and any possible new state's energy *H*~ν~ (using single-spin-flip dynamics) is 2*J* between the spin we choose to "flip" to move to the new state and that spin's neighbor. Thus, in a 1D Ising model, where each site has two neighbors (left and right), the maximum difference in energy would be 4 *J* .

## Running the Script!

Open terminal in your folder and execute the below command to create an executable file of your fortran code

```bash
gfortran Ising1D.f90 -o run.exe
```

Now, run the executable file

```bash
./run.exe
```

Cheers! Now, we proceed to plot our output using `gnuplot`. So, to do that, I have already prepared `.gnu` files, which contains the required `gnuplot` commands to get the desired output. Just execute the below commands in terminal to run them too. It will create two `.png` files of the desired plots!

```bash
gnuplot Plot-Mag.gnu
gnuplot Plot-Pot.gnu
```

Feel free to reach out with any questions or suggestions. I'm here to help!

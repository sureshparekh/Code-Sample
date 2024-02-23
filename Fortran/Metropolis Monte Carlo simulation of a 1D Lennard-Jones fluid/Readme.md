# Solving 1-D Lennard-Jones Potential using Metropolis Monte Carlo Method


In computational chemistry, molecular physics, and physical chemistry the **Lennard-Jones potential** is an intermolecular pair potential. Out of all the intermolecular potentials, the Lennard-Jones potential is probably the one that has been the most extensively studied. It is considered an archetype model for simple yet realistic intermolecular interactions.

The Lennard-Jones potential models soft repulsive and attractive (van der Waals) interactions. Hence, the Lennard-Jones potential describes electronically neutral atoms or molecules. The commonly used expression for the Lennard-Jones potential is

![{\displaystyle V_{\text{LJ}}(r)=4\varepsilon \left[\left({\frac {\sigma }{r}}\right)^{12}-\left({\frac {\sigma }{r}}\right)^{6}\right],}](https://wikimedia.org/api/rest_v1/media/math/render/svg/82d16f77cae964a5c4c52fb89165dd5d596ee03f)

where **r** is the distance between two interacting particles, **ε** is the depth of the potential well (usually referred to as 'dispersion energy'), and **σ** is the distance at which the particle-particle potential energy **V** is zero (often referred to as 'size of the particle'). The Lennard-Jones potential has its minimum at a distance of, ![{\displaystyle r=r_{\rm {min}}=2^{1/6}\sigma ,}](https://wikimedia.org/api/rest_v1/media/math/render/svg/2bb226b8a7c821063a18eaa90413f86eb90ae410)

 where the potential energy has the value ![{\displaystyle V=-\varepsilon .}](https://wikimedia.org/api/rest_v1/media/math/render/svg/ddd4f5f5e84d49724da899b3a3c3e8c9f85d67c0)

## Running the Script!

Open terminal in your folder and execute the below command to create an executable file of your fortran code

```bash
gfortran MCMC.f90 -o run.exe
```

Now, run the executable file

```bash
./run.exe
```

Cheers! Now, we proceed to plot our output using `gnuplot`. So, to do that, I have already prepared a `.gnu` file, which contains the required `gnuplot` commands to get the desired output. Just execute the below commands in terminal. It will create a `.png` files of the desired plot!

```bash
gnuplot Plot_graph.gnu
```

Feel free to reach out with any questions or suggestions. I'm here to help!

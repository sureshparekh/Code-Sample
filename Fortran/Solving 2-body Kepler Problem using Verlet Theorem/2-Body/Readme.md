# Solving 2-body Kepler's Problem using Verlet Theorem

In classical mechanics, the **Kepler problem** is a special case of the two-body problem, in which the two bodies interact by a central force *F* that varies in strength as the inverse square of the distance *r* between them. The force may be either attractive or repulsive. The problem is to find the position or speed of the two bodies over time given their masses, positions, and velocities. Using classical mechanics, the solution can be expressed as a Kepler Orbit using six orbital elements.

**Verlet integration Theorem** is a numerical method used to integrate Newton's equations of motion. It is frequently used to calculate trajectories of particles in molecular dynamics simulations and computer graphics.

To solve the problem, we begin with supplying the initial positions and velocities of the two particles to the code. Also, their masses. Furthermore, we proceed to calculate the instantaneous positions of both the particles alongwith their acceleration due to the mutual forces acting on the particles. We store the values of the positions and velocities of both the particles in seperate files. It can be used to plot the position and Velocity plots.

## Running the Script!

Open terminal in your folder and execute the below command to create an executable file of your fortran code

```bash
gfortran K2Body.f90 -o run.exe
```

Now, run the executable file

```bash
./run.exe
```

Cheers! We have generated the positions and velocities files of both the particles. Now, we proceed to plot them using `gnuplot`. So, to do that, I have already prepared `.gnu` files, which contains the required `gnuplot` commands to get the desired output. Just execute the below commands in terminal to run them too. It will create two `.png` files of the desired plots!

```bash
gnuplot Plot-X.gnu
gnuplot Plot-V.gnu
```

Feel free to reach out with any questions or suggestions. I'm here to help!

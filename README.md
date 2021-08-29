# Asteroid-Simulation
Simulation of Asteroids in a Sun and Jupiter system using FORTRAN

The code included involves :

1) potential.f90 -> For the Potential energy in a rotating frame in the presence of two massive bodies.
2) RK.f90, Euler.f90, Verlet.f90, ER.f90 -> Simulation of a single asteroid near L4 Lagrange Point with Runge-Kutta 4th order, Euler, Verlet and Euler-Richardson Algorithm and their respective Energy stability over time.
3) Traj.f90 -> The trajectory of a single asteroid in a rotating frame with Sun and Jupiter
4) Asteroids.f90 -> Simulation of the asteroid belt
5) AsteroidGenerate.ipynb -> Generating initial positions and velocities of all asteroids to be simulated.

An assumption of circular Jupiter orbit has been made throughout all simulations and the equations used are the solution to Restricted Three Body Problem.

Final Report : https://drive.google.com/file/d/1FgU0SFafr6RlhSSJONa3SGpSzjsURL1o/view?usp=sharing

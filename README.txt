Numerical Solution of the Cahn-Hilliard equation - Imran Marwat

The Cahn-Hilliard equation is a 2nd order differential equation that models the behaviour of phase separation in a binary mixture. A solution is found based on the finite difference method (see https://en.wikipedia.org/wiki/Finite_difference_method for more information). 

CH_Class.py contains all the methods for simulating and sampling a binary mixture (2D array) that evolves under the CH equation.
CH_Script.py contains instructions for running the methods above for various tasks.

To Run: python3 CH_Script.py <dimension> <initialvalue> <task>

-Using the finite difference method, an algorithm is created to find the value of the binary mixture, at a given point in the lattice, for the (n+1)th timestep based on the lattice values in the nth timestep. The solution in the time domain ultimately converges to a stable solution.

-Initial values of -0.5 or 0.5 lead to droplets of oil or water respectively forming in the second liquid. An initial value of 0 leads to a non-uniform mixture which eventually separates into a layer of oil and water each. 

-Available tasks are:

<task> = "viz": animates the solution of the equation. User needs to close animation window to end program.

<task> = “data”: simulates the solution for a large number of timesteps. Once the system reaches equilibrium the free energy is recorded at regular intervals. The result is plotted and shows that the system correctly simulates a physical system which acts to minimise its free energy. 


import numpy as np
import sys
import matplotlib.pylab as plt
from matplotlib.animation import FuncAnimation
from CH_Class import PDE

#Define time and space increments.
dt, dx = 2,1

#Create instance of PDE class and add random noise.
A = PDE(dt,dx,int(sys.argv[1]),float(sys.argv[2]))
A.Noise()

#If user requires visual simulation
if sys.argv[3]=='viz':
    
    #Update function
    def UpdatePlot(*args):
        image = ax.imshow(A.order_array)
        for i in range(50):
            A.Sweep()
        return image,
    
    #Create animation
    fig,ax = plt.subplots()
    image = ax.imshow(A.order_array)
    ani = FuncAnimation(fig,UpdatePlot,blit=True)
    plt.show()

#If user requires free energy data
elif sys.argv[3]=='data':
    
    #Lists for plotting
    FE_list = []
    time_list = []
    
    #Simulates the system for a number of iterations
    for i in range(30000):
        A.Sweep()
        if i>=500 and i%500==0:
            #Measures free energy at regular intervals
            FE_list.append(A.Free_Energy())
            time_list.append(i)
            print(i)

    #Plots free energy vs timestep
    plt.plot(time_list,FE_list)
    plt.xlabel("Timestep")
    plt.ylabel("Free Energy")
    plt.title("Plot of the Free Energy of the order parameter lattice vs Timestep")
    plt.show()


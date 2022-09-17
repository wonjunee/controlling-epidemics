import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import animation
import numpy as np
import csv
import sys

import matplotlib.patches as patches

def get_circle(center, radius):
    angle = np.linspace( 0 , 2 * np.pi , 150 )     
    x = radius * n1 * np.cos( angle ) + center[0] * n1 - 0.5
    y = radius * n2 * np.sin( angle ) + center[1] * n2 - 0.5
    return x,y

#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

def open_csv(filename,nt,n1,n2):
    A = np.fromfile(filename, dtype=np.float64)
    return A.reshape((nt,n2,n1))

directory="./data"

#--------------------------------------------------
#   Getting n1 and nt
#--------------------------------------------------

with open("{}/parameters.csv".format(directory)) as F:
    csvReader = csv.reader(F)
    for i in csvReader:
        n1 = int(i[0])
        n2 = int(i[1])
        nt = int(i[2])
        c0 = float(i[3])
        c1 = float(i[4])
        c2 = float(i[5])
        beta = float(i[6])
        gamma = float(i[7])

rho0 = open_csv("{}/rho0.csv".format(directory),nt,n1,n2)
rho1 = open_csv("{}/rho1.csv".format(directory),nt,n1,n2)
rho2 = open_csv("{}/rho2.csv".format(directory),nt,n1,n2)

#--------------------------------------------------
#   Getting n1 and nt
#--------------------------------------------------


fig, ax = plt.subplots(1,1)
ax.plot(np.sum(np.sum(rho0,axis=1),axis=1)/(n1*n2),label="S with control")
ax.plot(np.sum(np.sum(rho1,axis=1),axis=1)/(n1*n2),label="I with control")
ax.plot(np.sum(np.sum(rho2,axis=1),axis=1)/(n1*n2),label="R with control")
ax.plot(np.sum(np.sum(rho0+rho1+rho2,axis=1),axis=1)/(n1*n2),label="Total with control")
ax.set_xlabel("Time $t$")
ax.set_ylabel("The total number of people")
ax.legend()

plt.show()

#--------------------------------------------------
#   Create animation
#--------------------------------------------------

directory="./data"
rho0 = open_csv("{}/rho0.csv".format(directory),nt,n1,n2)
rho1 = open_csv("{}/rho1.csv".format(directory),nt,n1,n2)
rho2 = open_csv("{}/rho2.csv".format(directory),nt,n1,n2)

def get_rect(x,y,w,h):
    return x*n1-0.5-w*n1, y*n2-0.5-h*n2, w*n1*2, h*n2*2
def save_animation():
    # First set up the figure, the axis, and the plot element we want to animate
    num = 3
    
    fig, ax = plt.subplots(1,num,figsize=(9,4))
    

    cax0 = ax[0].imshow(rho0[0], cmap='inferno', origin='lower')
    cax1 = ax[1].imshow(rho1[0], cmap='inferno', origin='lower')
    cax2 = ax[2].imshow(rho2[0], cmap='inferno', origin='lower')


    ax[0].set_axis_off()
    ax[1].set_axis_off()
    ax[2].set_axis_off()
    plt.tight_layout()

    vmax0 = np.max(rho0)
    vmax1 = np.max(rho1)
    vmax2 = np.max(rho2)
    vmax  = max(vmax0, vmax1, vmax2)
    


    # animation function.  This is called sequentially
    def animate(n):

        fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, hspace=0, wspace=0.1)

        cax0.set_array(rho0[np.maximum(0,n-1)])
        cax1.set_array(rho1[np.maximum(0,n-1)])
        cax2.set_array(rho2[np.maximum(0,n-1)])


        rho0sum = np.sum(rho0[np.maximum(0,n-1)])/(n1*n2)
        rho1sum = np.sum(rho1[np.maximum(0,n-1)])/(n1*n2)
        rho2sum = np.sum(rho2[np.maximum(0,n-1)])/(n1*n2)

        # cax0.set_clim(0, vmax)
        cax0.set_clim(0, vmax0)
        cax1.set_clim(0, vmax1)
        cax2.set_clim(0, vmax2*3)

        ax[0].set_title("Susceptible: {:.4f}".format(rho0sum))
        ax[1].set_title("Infected: {:.4f}".format(rho1sum))
        ax[2].set_title("Recovered: {:.4f}".format(rho2sum))

        # cax0.set_clim(0, 10)
        plt.suptitle("$\\beta$={:.2f}, $\\gamma$={:.2}, t={:.2f}".format(beta,gamma,n/nt))

        # if(n==0):
        #     plt.show()
        return cax0, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt+1, interval=100, blit=True)
    anim.save("video.mp4", fps=10)

save_animation()
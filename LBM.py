import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm

#Definition of lattice constants
RHO=1 #Lattice density
ULB=0.04 #Lattice velocity
NJ=130 #Number of lattices in x-direction
NI=300 #Number of lattices in y-direction
NSTEP=40000 #Total number of time-steps
RE=100 #Reynolds Number
WI=25 #Write Interval - Number of time-steps after which plots and files are prepared

#Define obstacle in line 91
#Obstacle parameters
CX=NI // 4 #X-coordinate of Center of obstacle 
CY=NJ // 2 #Y-coordinate of Center of obstacle 
R=5 #Radius of circle or Edge-length of Square

#Major and minor and minor axes of Ellipse
B=5
A=1.5*B

#Airfoil
AIRFOIL=30 #Last two digits of airfoil (NACA00XX)
CHORD=30 #Chord length in terms of lattices
AOA=-5 #Angle of Attack (Positive->Anti-clockwise)

#Calculation of LBM parameters
NULB=(ULB*2*R)/RE #Viscosity 
OMEGA=1/(3.*NULB + 0.5) #Time-relaxation parameter




def macroscopic(fin,v): #Calculation of Macroscopic Quantities
    rho=np.sum(fin, axis=0) #Density
    p=rho/3  #Pressure
    u=np.zeros((2,NI,NJ)) 
    for i in range(9):
        u[0,:,:] += v[i,0] * fin[i, :, :]
        u[1,:,:] += v[i,1] * fin[i, :, :]
    u /= rho #Velocity
    return rho,u,p

def equilibrium(rho,u,v,t):   #Calculation of Equilibrium Function
    usqr = 3./2. * (u[0]**2 + u[1]**2)
    feq= np.zeros((9,NI,NJ))
    for k in range(9):
        vu = 3.0 * (v[k,0]*u[0,:,:] + v[k,1]*u[1,:,:])
        feq[k,:,:] = rho*t[k] * (1 + vu + 0.5*vu**2 - usqr)
    return feq

def inivel(d,x,y):  #Initial Velocity with slight perturbation
    return (1-d)*ULB* (1 + 1e-4* np.sin(y/(NJ-1)*2*np.pi))

def obstacle_circle(x,y):  
    return (x-CX)**2 + (y-CY)**2 < (R)**2  #Equation of Circle

def obstacle_square(x,y):
    return abs(x-CX+y-CY) + abs(x-CX-y+CY) < R  #Equation of Square

def obstacle_rotated_square(x,y):
    return abs(x-CX) + abs(y-CY) < (R/math.sqrt(2))  #Equation of Rotated Square

def obstacle_ellipse(x,y):
    return ((x-CX)/A)**2 +((y-CY)/B)**2 < 1  #Equation of Ellipse
    
def naca_AOA(x,y):
    AOA_RAD=(AOA*np.pi/180)
    y_u=(5*0.3*CHORD*(0.2969*(((x-CX)/CHORD)**0.5) - 0.126*((x-CX)/CHORD) - 0.3516*(((x-CX)/CHORD)**2) + 0.2843*(((x-CX)/CHORD)**3) - 0.1036*(((x-CX)/CHORD)**4)))
    y_l=-y_u
    Y_U=y_u*np.cos(AOA_RAD) - (x-CX)*np.sin(AOA_RAD)
    Y_L=y_l*np.cos(AOA_RAD) - (x-CX)*np.sin(AOA_RAD)
    return ((x-CX)>0) & ((x-CX)<(CHORD*np.cos(AOA_RAD))) & ((y-CY)>Y_L) & ((y-CY)<Y_U) #Equation of NACA00XX Airfoil


def LBM():
    v=np.array([[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]])
    t=np.array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

    col1 = np.array([0,1,2])  #Left wall of Lattice
    col2 = np.array([3,4,5])  #Interior of Lattice
    col3 = np.array([6,7,8])  #Right wall of Lattice
    
    vel = np.fromfunction(inivel, (2,NI,NJ))
    feq = equilibrium(RHO,vel,v,t)
    fin = feq.copy()
   
    obstacle = np.fromfunction(obstacle_circle, (NI,NJ))  #Change obstacle definition
    
    obstacle_fx=np.zeros(NSTEP)
    obstacle_fy=np.zeros(NSTEP)
    vorticity=np.zeros((NI,NJ))

        
    for l in range(NSTEP):
        #print("RUNNING",l)
        obstacle_fx[l]=0
        obstacle_fy[l]=0
        
        fin[col3,-1,:] = fin[col3,-2,:] #Outflow Condition on Right Wall

        rho,u ,p= macroscopic(fin,v)
    
        #Inlet Velocity at Left Wall
        u[:,0,:] = vel[:,0,:]
        rho[0,:] = 1./(1.-u[0,0,:]) * (np.sum(fin[col2,0,:],axis=0) + 2*np.sum(fin[col3,0,:],axis=0))
        feq = equilibrium(rho,u,v,t)
        fin[col1,0,:] = feq[col1,0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]
        
        fout = fin - OMEGA * (fin - feq) #Collision Step
    
        #Bounce-back condition on obstacle
        for i in range(9):
            fout[i, obstacle]=fin[8-i, obstacle]

        #Streaming Step    
        for i in range(9):
            fin[i,:,:] = np.roll(np.roll(fout[i,:,:],v[i,0],axis=0),v[i,1],axis=1)
            
        #Calculating forces on the obstacle    
        for i in range(1,NI-1):
            for j in range(1,NJ-1):
                if (obstacle[i-1,j-1] or obstacle[i-1,j] or obstacle[i-1,j+1] or obstacle[i,j-1] or obstacle[i,j+1] or obstacle[i+1,j-1] or obstacle[i+1,j] or obstacle[i+1,j+1]):
                   for k in range(9):
                        obstacle_fx[l]=obstacle_fx[l] + v[k,0]*(fin[k,i,j]+fout[8-k,i,j])
                        obstacle_fy[l]=obstacle_fy[l] + v[k,1]*(fin[k,i,j]+fout[8-k,i,j])  

        #Calculation of vorticity
        for i in range(NI):
            for j in range(NJ):
                if i==0 or j==0 or i==(NI-1) or j==(NJ-1):
                    vorticity[i,j]=0
                else:
                    vorticity[i,j] = u[1,i+1,j] - u[1,i-1,j] - u[0,i,j+1] + u[0,i,j-1]
                
        
        #Post-processing Plots
        if l % WI== 0:
            plt.clf() 
            plt.axis('off')
            plt.imshow(np.sqrt(u[0]**2+u[1]**2).transpose(),cmap=cm.jet)
            plt.savefig("vel_"+str(l)+".png")
        
        #Post-processing Files
        if l % WI == 0:
            f1=open("time_"+str(l)+".dat","w+")
            f1.write('VARIABLES = "X", "Y", "V", "P", "VORTICITY" \n')
            f1.write("ZONE I="+str(NI)+", J="+str(NJ)+" F=POINT \n")
            f1.write("SOLUTIONTIME="+str(l)+"\n")
            for a in range(NJ):
                for b in range(NI):
                    f1.write(str(b)+" "+str(a)+" "+str(np.sqrt(u[0,b,a]**2+u[1,b,a]**2))+" "+str(p[b,a])+" "+str(vorticity[b,a])+"\n")
            f1.close()

    #Printing Drag and Lift Forces (Only at end of simulation)
    f2=open("drag_lift.dat","w+")
    f2.write('VARIABLES = "T", "DRAG", "LIFT" \n')
    for a in range(NSTEP):
        f2.write(str(a)+" "+str(obstacle_fx[a])+" "+str(obstacle_fy[a])+"\n")
    f2.close()

LBM()





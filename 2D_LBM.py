# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 09:05:50 2021

@author: Dev Hathi
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm

#Definition of lattice constants
RHOLB=1 #Lattice density
ULB=0.04 #Lattice velocity
NJ=130 #Number of lattices in x-direction
NI=300 #Number of lattices in y-direction
NSTEPS=10000 #Total number of time-steps
RE=100 #Reynolds Number
WI=250 #Write Interval - Number of time-steps after which plots and files are prepared
WI_DL=50
DL_COUNT=0

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

C=np.array([[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]])
W=np.array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])
COL1 = np.array([0,1,2])  #Left wall of Lattice
COL2 = np.array([3,4,5])  #Interior of Lattice
COL3 = np.array([6,7,8])  #Right wall of Lattice

FIN=np.zeros((9,NI,NJ))
FEQ=np.zeros((9,NI,NJ))
FOUT=np.zeros((9,NI,NJ))
RHO=np.zeros((NI,NJ))
U=np.zeros((2,NI,NJ))
P=np.zeros((NI,NJ))
VORTICITY=np.zeros((NI,NJ))
DRAG=np.zeros(NSTEPS+1)
LIFT=np.zeros(NSTEPS+1)

def CALC_QUANTITIES(): #Calculation of Macroscopic Quantities
    global U,P,RHO
    RHO=np.sum(FIN, axis=0) #Density
    P=RHO/3  #Pressure
    U[:,:,:]=0
    for k in range(9):
        U[0,:,:] += C[k,0] * FIN[k, :, :]
        U[1,:,:] += C[k,1] * FIN[k, :, :]
    U /= RHO #Velocity

def CALC_EQUILIBRIUM():   #Calculation of Equilibrium Function
    USQR = 3./2. * (U[0]**2 + U[1]**2)
    for k in range(9):
        VU = 3.0 * (C[k,0]*U[0,:,:] + C[k,1]*U[1,:,:])
        FEQ[k,:,:] = RHO*W[k] * (1 + VU + 0.5*VU**2 - USQR)

def CALC_INIVEL(D,X,Y):  #Initial Velocity with slight perturbation
    return (1-D)*ULB* (1 + 1e-4* np.sin(Y/(NJ-1)*2*np.pi))

def CALC_VORTICITY():
    for i in range(NI):
        for j in range(NJ):
            if i==0 or j==0 or i==(NI-1) or j==(NJ-1):
                VORTICITY[i,j]=0
            else:
                VORTICITY[i,j] = (U[1,i+1,j] - U[1,i-1,j] - U[0,i,j+1] + U[0,i,j-1])/2
                
def CALC_DRAG_LIFT():
    global DL_COUNT
    DRAG[DL_COUNT]=0
    LIFT[DL_COUNT]=0
    for i in range(1,NI-1):
        for j in range(1,NJ-1):
            if (OBSTACLE[i-1,j-1] or OBSTACLE[i-1,j] or OBSTACLE[i-1,j+1] or OBSTACLE[i,j-1] or OBSTACLE[i,j+1] or OBSTACLE[i+1,j-1] or OBSTACLE[i+1,j] or OBSTACLE[i+1,j+1]):
                for k in range(9):
                    DRAG[DL_COUNT]=DRAG[DL_COUNT] + C[k,0]*(FIN[k,i,j]+FOUT[8-k,i,j])
                    LIFT[DL_COUNT]=LIFT[DL_COUNT] + C[k,1]*(FIN[k,i,j]+FOUT[8-k,i,j]) 
    DL_COUNT +=1

def OBSTACLE_CIRCLE(X,Y):  
    return (X-CX)**2 + (Y-CY)**2 < (R)**2  #Equation of Circle

def OBSTACLE_SQUARE(X,Y):
    return abs(X-CX+Y-CY) + abs(X-CX-Y+CY) < R  #Equation of Square

def OBSTACLE_ROTATED_SQUARE(X,Y):
    return abs(X-CX) + abs(Y-CY) < (R/math.sqrt(2))  #Equation of Rotated Square

def OBSTACLE_ELLIPSE(X,Y):
    return ((X-CX)/A)**2 +((Y-CY)/B)**2 < 1  #Equation of Ellipse
    
def OBSTACLE_AIRFOIL(X,Y):
    AOA_RAD=(AOA*np.pi/180)
    y_u=(5*0.3*CHORD*(0.2969*(((X-CX)/CHORD)**0.5) - 0.126*((X-CX)/CHORD) - 0.3516*(((X-CX)/CHORD)**2) + 0.2843*(((X-CX)/CHORD)**3) - 0.1036*(((X-CX)/CHORD)**4)))
    y_l=-y_u
    Y_U=y_u*np.cos(AOA_RAD) - (X-CX)*np.sin(AOA_RAD)
    Y_L=y_l*np.cos(AOA_RAD) - (X-CX)*np.sin(AOA_RAD)
    return ((X-CX)>0) & ((X-CX)<(CHORD*np.cos(AOA_RAD))) & ((Y-CY)>Y_L) & ((Y-CY)<Y_U) #Equation of NACA00XX Airfoil


def BC_OUTLET():
    FIN[COL3,-1,:] = FIN[COL3,-2,:]

def BC_INLET():
    U[:,0,:] = VEL[:,0,:]
    RHO[0,:] = 1./(1.-U[0,0,:]) * (np.sum(FIN[COL2,0,:],axis=0) + 2*np.sum(FIN[COL3,0,:],axis=0))
    CALC_EQUILIBRIUM()
    FIN[COL1,0,:] = FEQ[COL1,0,:] + FIN[[8,7,6],0,:] - FEQ[[8,7,6],0,:]

def BC_WALL():
    for k in range(9):
        FOUT[k, OBSTACLE]=FIN[8-k, OBSTACLE]
        
def COLLISION():
    global FOUT
    FOUT = FIN - OMEGA*(FIN-FEQ)

def STREAMING():
    for k in range(9):
        FIN[k,:,:] = np.roll(np.roll(FOUT[k,:,:],C[k,0],axis=0),C[k,1],axis=1)
        
def WRITE_FILES():
    f1=open("time_"+str(t)+".dat","w+")
    f1.write('VARIABLES = "X", "Y", "V", "P", "VORTICITY" \n')
    f1.write("ZONE I="+str(NI)+", J="+str(NJ)+" F=POINT \n")
    f1.write("SOLUTIONTIME="+str(t)+"\n")
    for j in range(NJ):
        for i in range(NI):
            f1.write(str(i)+" "+str(j)+" "+str(np.sqrt(U[0,i,j]**2+U[1,i,j]**2))+" "+str(P[i,j])+" "+str(VORTICITY[i,j])+"\n")
    f1.close()

def SAVE_PLOTS():
    plt.clf() 
    plt.axis('off')
    plt.imshow(np.sqrt(U[0]**2+U[1]**2).transpose(),cmap=cm.jet)
    plt.savefig("vel_"+str(t)+".png")
    
def WRITE_DRAG_LIFT():
    f2=open("drag_lift.dat","w+")
    f2.write('VARIABLES = "T", "DRAG", "LIFT" \n')
    for n in range(DL_COUNT):
        f2.write(str(n*WI_DL)+" "+str(DRAG[n])+" "+str(LIFT[n])+"\n")
    f2.close()    
    
VEL = np.fromfunction(CALC_INIVEL, (2,NI,NJ))
OBSTACLE = np.fromfunction(OBSTACLE_CIRCLE, (NI,NJ))

RHO[:,:] = RHOLB
U[:,:,:]=VEL[:,:,:]
CALC_EQUILIBRIUM()
FIN=FEQ.copy()


for t in range(NSTEPS+1):
    BC_OUTLET()
    CALC_QUANTITIES()
    BC_INLET()
    COLLISION()
    BC_WALL()
    STREAMING()
    if t % WI_DL==0:
        CALC_DRAG_LIFT()
    if t % WI ==0:
        CALC_VORTICITY()
        WRITE_FILES()
        SAVE_PLOTS()

WRITE_DRAG_LIFT()




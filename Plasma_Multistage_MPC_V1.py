
######################################################################################
# IMPORT TOOLS
######################################################################################

from casadi import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import random
from datetime import datetime
startTime = datetime.now()
import scipy.io as sio
import scipy
import os
import sys
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.base import BaseEstimator
import time
from sympy.solvers import solve as sympysolve
from sympy import Symbol


# Define time stamp to distinguish saved files
TimeStamp = datetime.now().strftime('%Y-%m-%d-%Hh%Mm%Ss')

# Define the directory in which the file is saved (change this!)
directory = '/Users/adbonzanini/Box Sync/Berkeley/Research/Explicit MPC Paper/DNN-MPC-Plasma'

# Switch some functionalities on or off
gpSwitch = 1            # GP correction on/off
OFswitch = 0            # Offset-free approach on/off
TrainOnly = 0           # Carry out training only (when validating/testing)
saveSwitch = 0          # Save outputs on/off

######################################################################################
# LOAD DATA
######################################################################################
print('\nLoading data and training GP...')
# DEFINE LAG L
L = 1
# Import training data
df_train=pd.read_csv(directory+'/Supporting-Data-Files/TrainingDataSimulation.csv') 
# Import test data
df_test=pd.read_csv(directory+'/Supporting-Data-Files/TestDataSimulation.csv')

#-999 are recognized as a title, so they are not included in the dataset

######################################################################################
# FORMAT INTO TRAINING AND TEST DATA
######################################################################################

# EXTRACT RELEVANT DATA
Xtrain = df_train.iloc[:,1:-2]
Ytrain = df_train.iloc[:,-2:]
Xtest = df_test.iloc[:,1:-2]
Ytest = df_test.iloc[:,-2:]

Xtrain = np.array(Xtrain)
Ytrain = np.array(Ytrain)
Xtest = np.array(Xtest)
Ytest = np.array(Ytest)



######################################################################################
# GAUSSIAN PROCESS REGRESSION
######################################################################################

kernel=C(1, (1e-3, 1e+3))*RBF(1,(1e-3, 1e+3))
# kernel=C(1.07, constant_value_bounds='fixed')**2*RBF(0.588,  length_scale_bounds='fixed')
alphaNoise1 = 1.
# kernel2=C(1, constant_value_bounds='fixed')*RBF(0.1,length_scale_bounds='fixed')
alphaNoise2 = 1.

# GP for Temperature (gpT)
gpT = GaussianProcessRegressor(kernel=kernel, alpha=alphaNoise1, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=5)
GPRT = gpT.fit(Xtrain,Ytrain[:,0])
hyperparametersT = str(GPRT.kernel_)
hyp_sT = float(hyperparametersT.split('**')[0])
hyp_lT = float(hyperparametersT.split('scale=')[1].split(')')[0])
LLT = GPRT.L_  #Extract the cholesky decomposition of K

# GP for Intensity (gpI)
gpI = GaussianProcessRegressor(kernel=kernel, alpha=alphaNoise2, optimizer='fmin_l_bfgs_b', n_restarts_optimizer=5)
GPRI = gpI.fit(Xtrain,Ytrain[:,1])
hyperparametersI = str(GPRI.kernel_)
hyp_sI = float(hyperparametersI.split('**')[0])
hyp_lI = float(hyperparametersI.split('scale=')[1].split(')')[0])
LLI = GPRI.L_  #Extract the cholesky decomposition of K


# Make the prediction on the test data
Ypred, sigma = GPRT.predict(Xtest, return_std=True)
Ypred2, sigma2 = GPRI.predict(Xtest, return_std=True)
# Ipred, sigmaI = GPRI.predict(Xtest, return_std=True)#, return_cov=True)
Tpred = Ypred#Ypred[:,0]
Ipred = Ypred2#Ypred[:,1]
sigmaT = sigma
sigmaI = sigma2

maxSigma = np.array([max(sigmaT), max(sigmaI)]).reshape(-1,1)


print('Done! \n')
# Exit and show plots if only training is required (do not enter MPC loops)
if TrainOnly ==1:
    plt.figure() 
    plt.subplot(2,1,1)
    plt.plot(Tpred, 'r')
    plt.plot(Ytest[:,0], 'b')
    plt.plot(Tpred+3*sigmaT, 'k--')
    plt.plot(Tpred-3*sigmaT, 'k--')
    plt.legend(['Predicted', 'Test', '99% CI'])
    plt.ylabel('$\Delta T$')
    plt.subplot(2,1,2)
    plt.plot(Ipred, 'r')
    plt.plot(Ytest[:,1], 'b')
    plt.plot(Ipred+3*sigmaI, 'k--')
    plt.plot(Ipred-3*sigmaI, 'k--')
    plt.ylabel('$\Delta I$')
    plt.legend(['Predicted', 'Test', '99% CI'])
    plt.show()
    sys.exit()

sio.savemat(directory+'/Supporting-Data-Files/hyperparameters.mat', {'hyp_sT':hyp_sT, 'hyp_lT': hyp_lT, 'hyp_sI':hyp_sI, 'hyp_lI': hyp_lI, 'alphaNoise1':alphaNoise1, 'alphaNoise2':alphaNoise2})


###########################################################################################################
# DEFINE KERNEL "MANUALLY"
###########################################################################################################

def casadi_RBF(mat1, mat2, sigma, l):
    #Convert numpy arrays into DM objects
    # mat1 = DM(mat1)
    # mat2 = DM(mat2)
    # Norm of each row of each of the two matrices
    trnorms1 = vertcat(*[mtimes(mat1[j, :], mat1[j, :].T) for j in range(mat1.shape[0])])
    trnorms2 = vertcat(*[mtimes(mat2[j, :], mat2[j, :].T) for j in range(mat2.shape[0])])
    
    #Multiply by a vector of ones of appropriate dimension to create two [mat1.shape[0] x mat2.shape[0]] matrices
    k1 = mtimes(trnorms1, np.ones((mat2.shape[0], 1)).T)
    k2 = mtimes(np.ones((mat1.shape[0], 1)), trnorms2.T)
    
    #Now, having x1^2 and x2^2, we can calculate the kernel function
    k = k1 + k2 -2 * mtimes(mat1, mat2.T) # k <-- (x1-x2)^2 = x1^2 -2x1x2 +x2^2
    k *= - 1./(2 * l**2) #k <-- (1/(2l))(x1-x2)^2
    return (sigma**2)*exp(k)   # k <-- s^2 exp((1/(2l))*(x1-x2)^2)


Xtrain = np.matrix(Xtrain)
Xtest = np.matrix(Xtest)

'''
Ytrain = np.matrix(Ytrain)
KK = casadi_RBF(Xtrain, Xtest, hyp_sT, hyp_lT)
'''


###########################################################################################################
# MPC PARAMETERS
###########################################################################################################

#Define cost parameters
Q = 1*np.array([[10, 0], [0, 0.1]])
R = 1*np.array([[5, 0], [0, 2]])

RR = 0*np.array([[5, 0], [0,5]])    # Penalize abrupt changes in the inputs
QQ = 0*Q                            # Penalize abrupt changes in states
PN = 0*Q                            # Terminal cost weight

nx = Q.shape[0]
nu = R.shape[0]

Np = 3                              # Prediction horizon
N = 25                              # Simulation horizon
N_robust = 1                        # Robust horizon for multistage MPC

#Initial point
yi = np.array([[0.],[0.]])

# Set point(s) and time(s) at which the reference changes
ysp1 = np.array([[2.], [2.]])
tChange = 5
ysp2 = np.array([[5.], [2.]])
tChange2 = 10000
ysp3 = np.array([[0.], [0.]])


###########################################################################################################
# VARIABLES & LINEAR MODEL MATRICES
###########################################################################################################
print('Initialize MPC problem...')
#Declare the variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
u1 = MX.sym('u1')
u2 = MX.sym('u2')
x = vertcat(x1, x2)
u = vertcat(u1, u2)
ss = MX.sym('xssVar', 4)
wNoise = MX.sym('wNoise', 2)

# Model Parameters (metal: Tss = 38; Iss = 9.7; Qss = 3; Pss = 3)
modelMatM = sio.loadmat(directory+'/Supporting-Data-Files/MIMOmodelMetal')
Areal = modelMatM['A']
Breal = modelMatM['B']

# Plant-model mismatch (glass: Tss = 38; Iss = 14; Qss = 3; Pss = 3)
modelMatG = sio.loadmat(directory+'/Supporting-Data-Files/MIMOmodelGlass')
A = modelMatG['A']
B = modelMatG['B']
C = modelMatG['C']
offset = 0*np.array([[0], [0]])     # Choose whether to add an artificial offset (to debug "measurements" of consecutive differences)

# Noise
# np.random.seed(1)
wDist = 0*np.random.normal(0, 0.1, (nx,Np+L+1))
wkDist = 0*np.random.normal(0, 0.1, (nx,N+L+1))



###########################################################################################################
# OFFSET-FREE AND ESTIMATION MATRICES (IF NEEDED)
###########################################################################################################

# Offset-free tracking parameters
Bd = np.eye(2,2)
Cd = np.zeros((2,2))
Haug = np.eye(2,2)
Aaug = np.vstack([np.hstack([np.eye(2,2)-A, -B]),np.hstack([np.dot(Haug,C), np.zeros((2,2))])])

# Steady-state Kalman Gain
Adist = np.vstack([np.hstack([A, Bd]), np.hstack([np.zeros((2,2)), np.eye(2,2)])])
Bdist = np.vstack([B, np.zeros((2,2))])
Cdist = np.hstack([C, Cd])
Qnoise = np.diag([0.1, 0.1, 0.1, 0.1])
Rnoise = np.diag([0, 0])
Pinf = scipy.linalg.solve_discrete_are(Adist.T, Cdist.T, Qnoise, Rnoise)
Kinf = np.dot(np.dot(Pinf, Cdist.T),np.linalg.inv(np.dot(np.dot(Cdist, Pinf), Cdist.T)+Rnoise))
print(Kinf)
# sys.exit()

#LQR gain
Plqr = scipy.linalg.solve_discrete_are(A, B, Q, R)
Klqr = -np.linalg.solve(R + np.dot(np.dot(B.T, Plqr), B), np.dot(np.dot(B.T, Plqr), A))

# Discrete Lyapunov Equation
Plyap = scipy.linalg.solve_discrete_lyapunov(A+np.dot(B, Klqr), Q+np.dot(np.dot(Klqr.T, R), Klqr))

Psi = R + np.dot(np.dot(B.T, Plyap), B)



###########################################################################################################
# TIGHTENED CONSTRAINTS
###########################################################################################################

# Load the constraints
if gpSwitch == 1:
    Constraints = sio.loadmat(directory+'/supporting-data-files/stochConstraints')
else:
    Constraints = sio.loadmat(directory+'/supporting-data-files/robustConstraintsWorstCase')

Xtight = Constraints['Xtight']
# Xtight[1][0][1,:] = [1, 0, 3.5]  # Override terminal constraint manually

# Assign the loaded constraints to variables
Ycon = Xtight[0][0]
Ucon = Constraints['Ucon']
# Ucon = np.array([[1, 0, 7], [-1, 0, 1.5], [0, 1, 2], [0, -1, 2.5]])   # Override input constraints manually

# Terminal constraint
Xf = Xtight[-1][0]



###########################################################################################################
# DEFINE CASADI FUNCTIONS
###########################################################################################################

# Dynamics and stage cost
xNext = mtimes(A,x)+mtimes(B,u) + wNoise
y = mtimes(C,x) 
Lstage = mtimes(mtimes((y-ss[0:2,0]).T,Q),(y-ss[0:2,0])) + mtimes(mtimes((u-ss[2:,0]).T,R),(u-ss[2:,0])) 

# Functions for nominal vs. real dynamics
F = Function('F', [x, u, wNoise, ss], [xNext, Lstage],['x','u', 'wNoise', 'ss'],['xNext', 'Lstage'])
Freal = Function('Freal', [x, u, wNoise, ss], [mtimes(Areal,x)+mtimes(Breal,u)+ wNoise+np.linalg.solve(C,offset), Lstage],['x','u', 'wNoise', 'ss'],['xNext', 'Lstage'])

# Variable inputs for test data (since these will change within the loop)
xTest = MX.sym('xTest', 1, L*4)


# CadADi function for temperature GP prediction
KsT = casadi_RBF(Xtrain, xTest, hyp_sT, hyp_lT)
LkT = solve(LLT, KsT)
Yout = mtimes(LkT.T, solve(LLT, Ytrain[:,0]))

KssT = casadi_RBF(xTest, xTest,hyp_sT, hyp_lT)
KinvT = inv(mtimes(LkT, LkT.T))
VarT = KssT - mtimes([KsT.T, KinvT, KsT])

Fpred = Function('Fpred', [xTest], [Yout],['xTest'],['Yout'])

# CadADi function for intensity GP prediction
KsI = casadi_RBF(Xtrain, xTest, hyp_sI, hyp_lI)
LkI = solve(LLI, KsI)
YoutI = mtimes(LkI.T, solve(LLI, Ytrain[:,1]))

KssI = casadi_RBF(xTest, xTest,hyp_sI, hyp_lI)
KinvI = inv(mtimes(LkI, LkI.T))
VarI = KssT - mtimes([KsI.T, KinvI, KsI])

FpredI = Function('FpredI', [xTest], [YoutI],['xTest'],['YoutI'])

# CadADi function for run-up data (we need to keep a history of L points to make the prediction at the next step)
DxReal = Freal(x=x, u=u, wNoise = wNoise, ss=[0,0,0,0])['xNext'] - F(x=x, u=u, wNoise = [0,0], ss=[0,0,0,0])['xNext']
FDx = Function('FDx', [x, u, wNoise], [DxReal], ['x', 'u', 'wNoise'], ['DxReal'])


print('Done! \n')
# Function for setpoint calculator (if offset-free operation is required)
def spCalculator(yspIn, spYc, Aaug, Bd, Cd, Ycon, Ucon): 
    print('-------------------------------------------------')
    print('Calculating Setpoint...')
    print('-------------------------------------------------')
    Haug = np.eye(2,2)
        
    # Create empty nlp
    spVar =[]
    lbVar =[]
    ubVar =[]
    Var0 =[]
    gg=[]
    lbgg=[]
    ubgg=[]
    
    # Create variables
    Xsp = MX.sym('Xsp', 2)
    spVar +=[Xsp]
    # lbVar +=[-Ycon[1,2], -Ycon[3,2]]
    lbVar +=  [-inf, -inf]
    # ubVar +=[Ycon[0,2], Ycon[2,2]]
    ubVar +=[inf, inf]
    Var0 +=[0, 0]
    
    Usp = MX.sym('Usp', 2)
    spVar +=[Usp]    
    lbVar +=[-inf, -inf]
    ubVar +=[inf, inf]
    # lbVar +=[-Ucon[1,2], -Ucon[3,2]]
    # ubVar +=[Ucon[0,2], Ucon[2,2]]
    Var0 +=[0, 0]
    
    Ysp = MX.sym('Ysp', 2)
    spVar +=[Ysp]
    lbVar +=[-inf, -inf]
    ubVar +=[inf, inf]
    Var0 +=[0, 0]
    
    gg+= [Ysp-mtimes(C,Xsp)]
    lbgg+=[0, 0]
    ubgg+=[0, 0]
        
    Baug = vertcat(mtimes(Bd, spYc), yspIn-mtimes(mtimes(Haug,Cd), spYc))
    
    gg+= [mtimes(Aaug,vertcat(Xsp,Usp))-Baug]
    lbgg+=[0, 0, 0, 0]
    ubgg+=[0, 0, 0, 0]

    # Create an NLP solver
    Qs = np.array([[1, 0], [0,1]]) 
    spProb = {'f': mtimes(mtimes((Ysp-vertcat(float(yspIn[0]), float(yspIn[1]))).T, Qs),(Ysp-vertcat(float(yspIn[0]), float(yspIn[1]))))  , 'x': vertcat(*spVar), 'g': vertcat(*gg)}
    sol_opts = {'ipopt.print_level':0, 'ipopt.max_cpu_time':1.5}
    solver = nlpsol('solver', 'ipopt', spProb, sol_opts)
    
    # Solve the NLP
    spSol = solver(x0=Var0, lbx= lbVar, ubx= ubVar, lbg=lbgg, ubg=ubgg)
    sp_opt = spSol['x'].full().flatten()
    sp_opt = np.array([[sp_opt[4]], [sp_opt[5]], [sp_opt[2]], [sp_opt[3]]])
    
    # Print solver/feasibility status 
    feasMeasure = solver.stats()['return_status']
    print('-------------------------------------------------')
    print(feasMeasure)
    print('-------------------------------------------------')
    
    # Baug = np.vstack([np.dot(Bd, spYc), yspIn-np.dot(np.dot(Haug,Cd), spYc)])
    # spCon = np.linalg.solve(Aaug, Baug)
    # ySP = np.dot(C, spCon[0:2])
    
    # spCon[0:2] = ySP
    # xxsp = np.linalg.solve(C, yspIn-spYc)
    # uu = np.linalg.solve(B, xxsp - np.dot(A, xxsp))
    # spCon = np.array([[float(yspIn[0])], [float(yspIn[1])], [float(uu[0])], [float(uu[1])]])
    
    return sp_opt
       

spYc = offset

sp = spCalculator(ysp1, spYc, Aaug, Bd, Cd, Ycon, Ucon)
yss = sp[0:2]
uss = sp[2:4]

yss = [float(yss[0]), float(yss[1])]
uss = [float(uss[0]), float(uss[1])]



###########################################################################################################
# INITIALIZE SCENARIOS FOR MULTISTAGE MPC
###########################################################################################################
scenario_idx = [0, 1]
N_scenarios = len(scenario_idx)
w_i = [1./N_scenarios]*N_scenarios

# Initialize vectors to store trajectories
y1S = np.zeros([N_scenarios, Np+1])
y2S = np.zeros([N_scenarios, Np+1])
u1S = np.zeros([N_scenarios, Np])
u2S = np.zeros([N_scenarios, Np])
gp1S = np.zeros([N_scenarios, Np])
gp2S = np.zeros([N_scenarios, Np])

print('Begin loop...')


###########################################################################################################
# RUN MODEL FOR L STEPS BEFORE STARTING THE MPC
###########################################################################################################
print '    Calculate Run-up Data...'

#Initialize
xRunUp = np.matrix(np.zeros((nx,N+L+1)))
xRunUpM = np.matrix(np.zeros((nx,N+L+1)))
DxRunUp = np.matrix(np.zeros((nx,N+L)))
uRunUp = np.matrix(np.zeros((nu,N+L)))
XRunUp = np.zeros((N+1, 4*L))
Ycorrection = np.zeros((nx,N+1))


xRunUp[:,0] = np.linalg.solve(C, yi)
xRunUpM[:,0] = np.linalg.solve(C, yi)

# Run model for L steps and store the difference b/w "measurements" and model

for j in range(0,L):
    uRunUp[:,j] = 0.001*np.ones((2,1))
    xRunUp[:, j+1] = Freal(x=xRunUp[:,j], u=uRunUp[:,j], wNoise = wDist[:,j], ss=[0,0, 0, 0])['xNext']
    xRunUpM[:, j+1] = F(x=xRunUpM[:,j], u=uRunUp[:,j], wNoise = np.zeros((2,1)), ss = [0, 0, 0, 0])['xNext']
    DxRunUp[:,j] = xRunUp[:, j+1]- xRunUpM[:, j+1]  #Essentially we are capturing the noise. In reality we would capture the deviation.
    # xRunUpM[:, j+1] = xRunUp[:, j+1] # Update model prediction

# Store run-up data in appropriate form --> L deviations and L inputs in the first row.
XRunUp[0, :] = vertcat(np.dot(C,DxRunUp[:, 0:L]), uRunUp[:,0:L]).T

#Run CasADi Function
prediction = [Fpred(XRunUp[0,:]), FpredI(XRunUp[0,:])]
Ycorrection[0][0] = prediction[0]
Ycorrection[1][0] = prediction[1]

print('    Done!')


#Pre-allocation for speed
uopt = np.zeros((nx, Np))
uOptSeq = np.zeros((nu,N))
fopt = np.zeros((N,1))
yTr = np.zeros((nx, N+1))
yTrPred = np.matrix(np.zeros((nx, N+1)))
yModel = np.matrix(np.zeros((nx, N+1)))

# Initialization
dhat = np.array(prediction).T
ssPlot=[[float(ysp1[0]), float(ysp1[1])]]
ssPlot += ssPlot
YcMat = []
xki = xRunUp[:,L]

#Convert to list
uRunUp=np.ndarray.tolist(uRunUp)
xRunUp=np.ndarray.tolist(xRunUp)
xRunUpM=np.ndarray.tolist(xRunUpM)
DxRunUp=np.ndarray.tolist(DxRunUp)
Ycorrection=np.ndarray.tolist(Ycorrection)


###########################################################################################################
# MPC LOOP
###########################################################################################################
Tstart = time.time()
# MPC LOOP
for k in range(0, N):
    # if k==1:
    #     sys.exit()
    xhati = xki
    yki = np.dot(C,xki)
    yhati = yki
    yTr[:, 0]=np.transpose(yki)
    yTrPred[:,0] = yki
    Jactual = 0;
    gp1_opt = np.zeros((1, Np+1))
    gp2_opt = np.zeros((1, Np+1))
    
    
    # Update the initial condition of the MPC run-up arrays
    uRunUpMPC = [vertcat(uRunUp[0][k:L+k], uRunUp[1][k:L+k])]
    # Complete model predictions (nominal)
    xRunUpMPC = [vertcat(xRunUp[0][k:k+L+1], xRunUp[1][k:k+L+1])]
    xRunUpMPC = [[float(xRunUpMPC[0][0])]+[float(xRunUpMPC[0][2])], [float(xRunUpMPC[0][1])]+[float(xRunUpMPC[0][3])]]
    # Nominal model predictions (nominal)
    xRunUpM_MPC = [vertcat(xRunUpM[0][k:k+L+1], xRunUpM[1][k:k+L+1])]
    xRunUpM_MPC = [[float(xRunUpM_MPC[0][0])]+[float(xRunUpM_MPC[0][2])], [float(xRunUpM_MPC[0][1])]+[float(xRunUpM_MPC[0][3])]]
    DxRunUpMPC = [vertcat(DxRunUp[0][k:k+L], DxRunUp[1][k:k+L])]
    YcorrectionMPC = [vertcat(Ycorrection[0][k], Ycorrection[1][k])]


    print '\n\n################################# NEW OPTIMIZATION #################################\n\n'
    
    #   At each step k the entire OCP for all scenarios is solved!
    #   Therefore, we need to intialize empty variables for each
    #   step k.
    
    # Start with an empty NLP
    w=[]    #Array of all the variables we will be optimizing over
    w0 = []
    lbw = []
    ubw = []
    J = 0
    g=[]
    lbg = []
    ubg = []
    
    print(yki)
    
    
    #     "Lift" initial conditions. Note that the initial node
    #     is the same for all scenarios, so the double index is not
    #     required.
    
    Xk = MX.sym('X0', nx)
    w += [Xk]
    lbw += [float(xhati[0]), float(xhati[1])]
    ubw += [float(xhati[0]), float(xhati[1])]
    w0 += [float(xhati[0]), float(xhati[1])]
    
    Yk = MX.sym('Y0', nx)
    w += [Yk]
    lbw += [-inf, -inf]
    ubw += [inf, inf]
    w0 += [0, 0]      


    ###########################################################################################################
    # MPC LOOP FOR DIFFERENT SCENARIOS
    ###########################################################################################################
    for n_sc in range(0, N_scenarios):
        
        np.random.seed(n_sc)
        # wReal = 0*vertcat(np.random.uniform(-3.0*0.25, 3.0*0.25, (1,N+1)),np.random.uniform(-3.0*0.25, 3*0.25, (1, N+1)))
        wReal = vertcat(np.random.normal(0., 0.6, (1,N+1)),np.random.normal(0., 0.6, (1, N+1)))
    

        ###########################################################################################################
        # OPTIMAL CONTROL PROBLEM - OPEN LOOP OPTIMIZATION
        ###########################################################################################################
        
        # Formulate the NLP
        for i in range(Np):
            # New NLP variable for the control
            Uk = MX.sym('U_' + str(i)+'_'+str(n_sc), nu)
            w   += [Uk]
            lbw += [Ucon[1,2]/Ucon[1,0], Ucon[3,2]/Ucon[3,1]]
            ubw += [Ucon[0,2]/Ucon[0,0], Ucon[2,2]/Ucon[2,1]]
            w0  += [0]*nu
            
            
            YGP = MX.sym('YGP_' + str(i) +'_'+ str(n_sc), nx)
            w   += [YGP]
            lbw += [-9999*gpSwitch, -9999*gpSwitch]
            ubw += [9999*gpSwitch, 9999*gpSwitch]
            w0  += [0]*nx
            
            
            # Integrate until the end of the interval
            Fk = F(x=Xk, u=Uk, wNoise = gpSwitch*(YGP+scenario_idx[n_sc]*3*maxSigma), ss=yss+uss)
            Xk_end = Fk['xNext']
            # Yk_end = mtimes(C, Xk_end)+0*YGP
            J=J+w_i[n_sc]*Fk['Lstage']
            # Penalize abrupt changes
            J = J + mtimes(mtimes((Uk-uopt[:,i]).T, RR), Uk-uopt[:,i]) #+ mtimes(mtimes((Yk_end-Yk).T, QQ), Yk_end-Yk)


            #Predict deviation from model using GP and the L last steps: x(k+1) = Ax + Bu + g(x,u)       
            uRunUpMPC += [Uk]
            uRunUpMPC = uRunUpMPC[-L:]
            
            xRunUpMPC = xRunUpMPC[-L:]
            xRunUpM_MPC = xRunUpM_MPC[-L:]

            xRunUpMPC += [F(xRunUpMPC[0], uRunUpMPC[0], wDist[:,i], [0,0,0,0])[0]]#<======== !!
            xRunUpM_MPC += [F(xRunUpMPC[0], uRunUpMPC[0], np.zeros((2,1)), [0,0,0,0])[0]]
            
            DxRunUpMPC += [mtimes(C, xRunUpMPC[L] - xRunUpM_MPC[L]) +YcorrectionMPC[i]]  #<==== TECHNICALLY Dy not Dx!!  
            DxRunUpMPC = DxRunUpMPC[-L:]
            
            xx = DxRunUpMPC+uRunUpMPC
            xx = reshape(vertcat(*xx), 1, -1)

            #Run CasADi Function
            if i==0:
                YcorrectionMPC += [vertcat(Fpred(xx), FpredI(xx))] # Correct for first prediction
            else:
                YcorrectionMPC += [0*YcorrectionMPC[1]] # Correct for future predictions
                #YcorrectionMPC += [vertcat(Fpred(xx), FpredI(xx))]
            
            g   += [Yk-mtimes(C,Xk)]
            lbg += [0]*nx
            ubg += [0]*nx
            
            g += [gpSwitch*YGP-gpSwitch*YcorrectionMPC[i+1]] 
            lbg += [0]*nx
            ubg += [0]*nx 
 
            # New NLP variable for state at end of interval
            Xk = MX.sym('X_' + str(i+1)+'_'+ str(n_sc), nx)
            w   += [Xk]
            lbw += [-inf, -inf]
            ubw += [inf, inf]
            w0  += [0]*nx
            
            Yk = MX.sym('Y_' + str(i+1)+'_'+ str(n_sc), nx)
            # Ycon = Xtight[i+1][0]             #<------------- TIGHT CONSTRAINTS
            Ycon = Xtight[0][0][[0, 2,1,3],:] #<------------- ORIGINIAL CONSTRAINTS
            w   += [Yk]
            lbw += [-inf, -inf]
            ubw += [inf, inf]
            w0  += [0]*nx
            
            # Add equality constraints
            g   += [mtimes(Ycon[:,0:2], Yk)]
            lbg += [-inf]*len(Ycon[:,2])
            ubg += np.ndarray.tolist(Ycon[:,2])
            
            g   += [Xk_end-Xk]
            lbg += [0]*nx
            ubg += [0]*nx
                

        # Terminal cost and constraints (Xk --> i+1)
        # Terminal Cost
        J = J + w_i[n_sc]*mtimes(mtimes((Yk-vertcat(yss)).T,PN),(Yk-vertcat(yss)))

        # # Terminal Constraint
        g += [mtimes(Xf[:,0:2], Yk-yss)]
        # g+=[Yk]
        lbg += [-inf]*len(Xf[:,2])
        ubg += np.ndarray.tolist(Xf[:,2])
        
        # Equality constraint to make sure that Yk at the last step is equal to  C*Xk
        g += [Yk-mtimes(C,Xk)]
        lbg += [0]*nx
        ubg += [0]*nx
        
    # Non-anticipativity constraints
    Nrepeat = (len(w)-2)/N_scenarios
    NstepVars = 4
    for con_idx1 in range(2, len(w)-2-Nrepeat, 4):
        for con_idx2 in range(con_idx1, len(w)-2-Nrepeat, Nrepeat):
            g += [w[con_idx2]-w[con_idx2+Nrepeat]]
            lbg += [0]*nx
            ubg += [0]*nx        
        
    # Create an NLP solver
    prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
    # sol_opts ={}
    sol_opts = {'ipopt.print_level':0, 'ipopt.max_cpu_time':10.}
    solver = nlpsol('solver', 'ipopt', prob, sol_opts)
    # solver = qpsol('solver', 'qpoases', prob)
        
    # Solve the NLP
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    w_opt = sol['x'].full().flatten()
    J_opt = sol['f'].full().flatten()
    
    y1_opt = w_opt[2::8]
    y2_opt = w_opt[3::8]
    u1_opt = w_opt[4::8]
    u2_opt = w_opt[5::8]
    gp1_opt = w_opt[6::8]
    gp2_opt = w_opt[7::8]
     
    feasMeasure = solver.stats()['return_status']
    print('----------------------------------------------')
    print(feasMeasure)
    print('----------------------------------------------')
    
    # Assign optimization outputs into the respective scenarios
    for ii in range(0, N_scenarios):
        y1S[ii,:] = np.hstack([y1_opt[0], np.split(y1_opt[1:], N_scenarios)[ii]])
        y2S[ii,:] = np.hstack([y2_opt[0], np.split(y2_opt[1:], N_scenarios)[ii]])
        u1S[ii,:] = np.split(u1_opt, N_scenarios)[ii]
        u2S[ii,:] = np.split(u2_opt, N_scenarios)[ii]
        gp1S[ii,:] = np.split(gp1_opt, N_scenarios)[ii]
        gp2S[ii,:] = np.split(gp2_opt, N_scenarios)[ii]
    
    
    # Perform plant simulations for the next step using the plant model
    uopt = np.vstack([u1_opt, u2_opt])
    uOptSeq[:,k] = uopt[:,0]
    fopt[k] = J_opt
    
    # If this is the first simulation step, define Ycorrection from solution 
    if k==0:
        YcorrSol = np.array([[gp1_opt[0]],[gp2_opt[0]]])
        
    # Otherwise, stack the solution to the existing matrix
    else:    
        YcorrSol = np.hstack([YcorrSol, np.array([[gp1_opt[0]],[gp2_opt[0]]])])

    
    #     Update intial condition. Can change the if statement to define 
    #     when the plant-model mismatch is introduced (e.g. glass-to-metal
    #     transition).
    if k<0:
        xki = mtimes(A,xki)+mtimes(B,uopt[:,0]) + wReal[:,k]
        yki = mtimes(C,xki) + offset 
    else:
        xki = mtimes(Areal,xki)+mtimes(Breal,uopt[:,0]) + wReal[:,k]
        yki = mtimes(C,xki) + offset 
    
    #     Disturbance estimator/Observer. Note that due to the state feedback
    #     assumption, xhat is not used, unless the offset-free formulation is
    #     implemented (see if statement below, where xhat is overwritten if 
    #     the GP formulation is implemented).
    x_d_hat = mtimes(Adist, vertcat(xhati, dhat)) +mtimes(Bdist, uopt[:,0]) + np.dot(Kinf, yki - np.dot(C, xhati)-np.dot(Cd, dhat))
    xhati = x_d_hat[0:2]
    dhat = x_d_hat[2:]
    
    if OFswitch == 1:
        Ycorrection[0][k+1] = float(dhat[0])
        Ycorrection[1][k+1] = float(dhat[1])
        
    elif gpSwitch == 1:
        # Disturbance "measurement" (comment out for disturbance estimation)
        # dhat = yki - (mtimes(A,xhati)+mtimes(B,uopt[:,0]))
        
        #     State estimate is given by applying the GP correction (estimation
        #     is overwritten.
        xhati = mtimes(A, xhati) +mtimes(B, uopt[:,0])+YcorrSol[:,k]
        yModel[:,k+1] = mtimes(C, xhati)

        #     Disturbance estimate used to initialize the plant-model mismatch 
        #     at the next simulation step
        Ycorrection[0][k+1] = float(dhat[0])
        Ycorrection[1][k+1] = float(dhat[1])
        
    else:
        Ycorrection[0][k+1] = 0
        Ycorrection[1][k+1] = 0
    
    
    # State Feedback (comment out for output feedback)!!
    xhati = yki 
        
    xPred = np.array((mtimes(A,xki) + mtimes(B,uopt[:,0])) + 0*np.vstack((Ycorrection[0][k], Ycorrection[1][k]))).reshape(2,)
    yTrPred[:,k+1] = mtimes(C, xPred)
    
    yTr[:,k+1] = np.array(yki).reshape(2,)
    # yTr[:,k+1] = np.array(xhati).reshape(2,)
    
    #Update "real" run-up data
    uRunUp[0][k+L] = uopt[0,0]
    uRunUp[1][k+L] = uopt[1,0]
    xRunUp[0][k+L+1] = float(xhati[0])
    xRunUp[1][k+L+1] = float(xhati[1]) 
    xRunUpM[0][k+L+1] = float(yTrPred[0,k+1])
    xRunUpM[1][k+L+1] = float(yTrPred[1,k+1])
    DxRunUp[0][k+L] = float(FDx(xhati, uopt[:,0], YcorrSol[:,-1]+wkDist[:,L+k])[0])
    DxRunUp[1][k+L] = float(FDx(xhati, uopt[:,0], YcorrSol[:,-1]+wkDist[:,L+k])[1])
    
    DyRunUp = np.ndarray.tolist(np.dot(C, np.matrix([DxRunUp[0][1], DxRunUp[1][1]]).T))
    xxReal = [DyRunUp[0][0], DyRunUp[1][0], u1_opt[0], u2_opt[0]] #<==

    # YcorrNext = [Fpred(xxReal), FpredI(xxReal)]
    # Ycorrection[0][k+1] = float(gp1_opt[0]) #<=== predicted =/= real
    # Ycorrection[1][k+1] = float(gp2_opt[0])
    
    # Ycorrection[0][k+1] = DyRunUp[0][0]
    # Ycorrection[1][k+1] = DyRunUp[1][0]
    
    YcMat +=[np.vstack([gp1_opt, gp2_opt])]
    
    # No offset-free formulation => no setpoint correction
    if OFswitch == 0:
        spYc = np.zeros((2,1))
        # spYc = 0*np.array([[gp1_opt[0]], [gp2_opt[0]]])  #<==== update sp at the beginning or at the end?
        
    # Otherwise => correct by the amount of the estimated disturbance
    else:
        spYc = np.array([[Ycorrection[0][k+1]], [Ycorrection[1][k+1]]])
        # spYc = np.array([[dhat[0]], [dhat[1]]])

    if k>=0 and k<tChange-1:
        # yss = np.ndarray.tolist(ysp1.T)[0]
        sp = spCalculator(ysp1, spYc, Aaug, Bd, Cd, Ycon, Ucon)
        yss = sp[0:2]
        uss = sp[2:4]
        yss = [float(yss[0]), float(yss[1])]
        uss = [float(uss[0]), float(uss[1])]
        # uss = spCalculator(yss, spYc)
        # uss = [float(uss[0]), float(uss[1])]   
        ssPlot+=[[float(ysp1[0]), float(ysp1[1])]]  
    elif k>=tChange-1 and k<tChange2-1:
        sp = spCalculator(ysp2, spYc, Aaug, Bd, Cd, Ycon, Ucon)
        yss = sp[0:2]
        uss = sp[2:4]
        yss = [float(yss[0]), float(yss[1])]
        uss = [float(uss[0]), float(uss[1])]
        ssPlot+=[[float(ysp2[0]), float(ysp2[1])]] 
    elif k>=tChange2-1:
        sp = spCalculator(ysp3, spYc, Aaug, Bd, Cd, Ycon, Ucon)
        yss = sp[0:2]
        uss = sp[2:4]
        yss = [float(yss[0]), float(yss[1])]
        uss = [float(uss[0]), float(uss[1])]
        ssPlot+=[[float(ysp3[0]), float(ysp3[1])]] 
    
    

###########################################################################################################       
Tend = time.time()
print '\n'
print 'time =', Tend-Tstart,'s'

'''
# Constraint violation fraction
CVF = yTr[0][tChange+20:]>=5
CVF = CVF.astype(int)
if len(CVF)>=1:
    CVF = float(sum(CVF))*100./float(len(CVF))
    CVFvec+=[CVF]

CVFavg = np.mean(CVFvec) 
print('Constraint Violation = ', CVFavg, '%')
'''
if saveSwitch ==1:
    sio.savemat('/users/adbonzanini/Box Sync/Berkeley/Research/Gaussian Process/MIMO GP State Feedback Substrates/Output Data Files/'+TimeStamp+'_trajectories.mat', {'xTraj':yTr, 'uTraj': uOptSeq, 'xsp':ssPlot, 'Q':Q,'R':R, 'RR':RR, 'Np':Np, 'N':N, 'gpSwitch':gpSwitch, 'OFswitch':OFswitch, 'Ucon':Ucon, 'Xtight':Xtight, 'y1MC':y1MC, 'y2MC':y2MC, 'u1MC':u1MC, 'u2MC':u2MC})


# Plot figures
plt.figure()
plt.subplot(2,1,1)
plt.plot(yTr[0,:], 'r-')
plt.step(np.linspace(0,N+1, N+2), np.array(ssPlot)[:,0], 'k-')
plt.plot([0,N], [Xtight[0][0][0,2]/Xtight[0][0][0,0]]*2, 'k--')
plt.plot([0,N], [Xtight[0][0][2,2]/Xtight[0][0][2,0]]*2, 'k--')
plt.ylabel("$T$")
plt.xlabel("$k$")
plt.ylim([-2, Xtight[0][0][0,2]+1.5])
plt.subplot(2, 1, 2)
plt.plot(yTr[1,:], 'r-')
plt.step(np.linspace(0,N+1, N+2), np.array(ssPlot)[:,1], 'k-')
plt.plot([0,N], [Xtight[0][0][1,2]/Xtight[0][0][1,1]]*2, 'k--')
plt.plot([0,N], [Xtight[0][0][3,2]/Xtight[0][0][3,1]]*2, 'k--')
plt.ylabel("$I$")
plt.xlabel("$k$")
plt.ylim([-1.5, Xtight[0][0][1,2]+0.5])

plt.figure()
plt.subplot(2,1,1)
plt.plot(uOptSeq[0,:], 'b')
plt.ylabel("$Flow (q)$")
plt.xlabel("$k$")
plt.ylim(np.min(uOptSeq[0,:])-1, np.max(uOptSeq[0,:])+1)
plt.plot([0,N], [Ucon[0,2]/Ucon[0,0]]*2, 'k--')
plt.plot([0,N], [Ucon[1,2]/Ucon[1,0]]*2, 'k--')
plt.ylim([-Ucon[1,2]-0.5, Ucon[0,2]+0.5])
plt.subplot(2, 1, 2)
plt.plot(uOptSeq[1,:], 'b')
plt.ylabel("$Power (P)$")
plt.xlabel("$k$")
plt.ylim(np.min(uOptSeq[1,:])-1, np.max(uOptSeq[1,:])+1)
plt.plot([0,N], [Ucon[2,2]/Ucon[2,1]]*2, 'k--')
plt.plot([0,N], [Ucon[3,2]/Ucon[3,1]]*2, 'k--')
plt.ylim([-Ucon[3,2]-0.5, Ucon[2,2]+0.5])
plt.show()





'''
plt.figure()  
plt.subplot(2,1,1) 
plt.plot(range(0,N+1), xMC[0,:], 'r-')
plt.plot(ssPlot, 'g-')
for k in range(1,N+1):
    plt.plot(range(0,N+1), xMC[1,:], 'r-')
plt.plot([0,N], [Xtight[0][0][0,1]/Xtight[0][0][0,0]]*2, 'k--')
plt.plot([0,N], [Xtight[0][0][1,1]/Xtight[0][0][1,0]]*2, 'k--')
plt.xlim([0, N+1])
plt.ylim([-6.5, 6.5])
plt.ylabel("$x_k$")
plt.xlabel("$k$")
plt.legend(['$x$', '$x_{sp}$',])# 'x_{predicted}'])
# plt.plot(range(0, N+1), yTrPred[0,:], 'g--')

# plt.subplot(2,1,2)
# plt.plot(range(0,N), uMC[0,:], 'b-')
# plt.plot(range(0,N), uMC[1,:], 'g-')
# plt.plot([0,N], [Ucon[0,1]/Ucon[0,0]]*2, 'k--')
# plt.plot([0,N], [Ucon[1,1]/Ucon[1,0]]*2, 'k--')
# plt.xlim([0, N+1])
# plt.ylim([-5, 5])
# plt.ylabel("$u_k$")
# plt.xlabel("$k$")

plt.show()

'''



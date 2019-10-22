#
#   This file is part of do-mpc
#
#   do-mpc: An environment for the easy, modular and efficient implementation of
#        robust nonlinear model predictive control
#
#   Copyright (c) 2014-2016 Sergio Lucia, Alexandru Tatulea-Codrean
#                        TU Dortmund. All rights reserved
#
#   do-mpc is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as
#   published by the Free Software Foundation, either version 3
#   of the License, or (at your option) any later version.
#
#   do-mpc is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with do-mpc.  If not, see <http://www.gnu.org/licenses/>.
#

# Additional imports
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


# This is the main path of your do-mpc installation relative to the execution folder
path_do_mpc = '/users/adbonzanini/Box Sync/Berkeley/Research/Explicit MPC Paper/do-mpc/'
print path_do_mpc
# Add do-mpc path to the current directory
import sys
sys.path.insert(0,path_do_mpc+'code')
# Do not write bytecode to maintain clean directories
sys.dont_write_bytecode = True
# Compatibility for python 2.7 and python 3.0
from builtins import input
# Start CasADi
from casadi import *
# Import do-mpc core functionalities
import core_do_mpc
# Import do-mpc plotting and data managament functions
import data_do_mpc


"""
-----------------------------------------------
Train the GP model
-----------------------------------------------
"""

######################################################################################
# IMPORT DATA AND DEFINE PARAMETERS
######################################################################################

# Import training data
df_train=pd.read_csv('/Users/adbonzanini/Box Sync/Berkeley/Research/Gaussian Process/MIMO GP State Feedback Substrates/Supporting Data Files/TrainingDataSimulation.csv') 
# Import test data
df_test=pd.read_csv('/Users/adbonzanini/Box Sync/Berkeley/Research/Gaussian Process/MIMO GP State Feedback Substrates/Supporting Data Files/TestDataSimulation.csv')
# Define how mant "history points" we want to be keeping track of
L = 1
# Initialize counter
j=0

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

#'''
#Plot results
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
#'''

print 'GP Trained!'
sio.savemat('/users/adbonzanini/Box Sync/Berkeley/Research/Explicit MPC Paper/DNN-MPC-Plasma/Supporting-Data-Files/hyperparameters.mat', {'hyp_sT':hyp_sT, 'hyp_lT': hyp_lT, 'hyp_sI':hyp_sI, 'hyp_lI': hyp_lI, 'alphaNoise1':alphaNoise1, 'alphaNoise2':alphaNoise2})
plt.show()


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



###########################################################################################################
# CASADI FUNCTIONS FOR GP PREDICTIONS
###########################################################################################################

# Variable test inputs instead of test data (test inputs are going to change at every simulation step)
xTest = MX.sym('xTest', 1, L*4)

# Prediction for temperature
KsT = casadi_RBF(Xtrain, xTest, hyp_sT, hyp_lT)
LkT = solve(LLT, KsT)
Yout = mtimes(LkT.T, solve(LLT, Ytrain[:,0]))
Fpred = Function('Fpred', [xTest], [Yout],['xTest'],['Yout'])

# Prediction for intensity
KsI = casadi_RBF(Xtrain, xTest, hyp_sI, hyp_lI)
LkI = solve(LLI, KsI)
YoutI = mtimes(LkI.T, solve(LLI, Ytrain[:,1]))
FpredI = Function('FpredI', [xTest], [YoutI],['xTest'],['YoutI'])

"""
-----------------------------------------------
do-mpc: Definition of the do-mpc configuration
-----------------------------------------------
"""

# Import the user defined modules
import template_model
import template_optimizer
import template_observer
import template_simulator


# Create the objects for each module
model_1 = template_model.model()
# Create an optimizer object based on the template and a model
optimizer_1 = template_optimizer.optimizer(model_1)
# Create an observer object based on the template and a model
observer_1 = template_observer.observer(model_1)
# Create a simulator object based on the template and a model
simulator_1 = template_simulator.simulator(model_1)
# Create a configuration
configuration_1 = core_do_mpc.configuration(model_1, optimizer_1, observer_1, simulator_1)

# Set up the solvers
configuration_1.setup_solver()



"""
----------------------------
do-mpc: MPC loop
----------------------------
"""
while (configuration_1.simulator.t0_sim + configuration_1.simulator.t_step_simulator < configuration_1.optimizer.t_end):

    """
    ----------------------------
    do-mpc: Optimizer
    ----------------------------
    """
    # Make one optimizer step (solve the NLP)
    configuration_1.make_step_optimizer()

    """
    ----------------------------
    do-mpc: Simulator
    ----------------------------
    """
    # Simulate the system one step using the solution obtained in the optimization
    configuration_1.make_step_simulator()

    """
    ----------------------------
    do-mpc: Observer
    ----------------------------
    """
    # Make one observer step
    configuration_1.make_step_observer()

    """
    ------------------------------------------------------
    do-mpc: Prepare next iteration and store information
    ------------------------------------------------------
    """
    # Store the information
    configuration_1.store_mpc_data()

    # Set initial condition constraint for the next iteration
    configuration_1.prepare_next_iter()

    """
    ------------------------------------------------------
    do-mpc: Plot MPC animation if chosen by the user
    ------------------------------------------------------
    """
    # Plot animation if chosen in by the user
    data_do_mpc.plot_animation(configuration_1)
    
    # Keep the history of the L most recent points
    xopt = configuration_1.mpc_data.mpc_states
    uopt = configuration_1.mpc_data.mpc_control
    
    # Collect and concatenate past states
    xplusNominal = mtimes(A, xopt[-L-1:-L,:].T)+mtimes(B, uopt[-L-1:-L,:].T)
    xDiff = xopt[-L:,:].T - xplusNominal
    xHistory = vertcat(xDiff)
    # Collect and concatenate past inputs
    uHistory = vertcat(*uopt[-L:,:])
    xx = vertcat(xHistory, uHistory).T
    
    print(xplusNominal)
    print(xDiff)
    print(xHistory)
    
    # GP Predictions
    w1 = vertcat(Fpred(xx))
    w2 = vertcat(FpredI(xx))  
    print(w1)   
    print(w2)
    
    j = j+1
    def globalVariables():
        return j, L
    if j==3:
        sys.exit()

"""
------------------------------------------------------
do-mpc: Plot the closed-loop results
------------------------------------------------------
"""

data_do_mpc.plot_mpc(configuration_1)

# Export to matlab if wanted
data_do_mpc.export_to_matlab(configuration_1)


input("Press Enter to exit do-mpc...")

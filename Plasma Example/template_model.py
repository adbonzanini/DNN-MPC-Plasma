#-*- coding: utf-8 -*-
#
#   This file is part of do-mpc
#
#   do-mpc: An environment for the easy, modular and efficient implementation of
#        robust nonlinear model predictive control
#
#   Copyright (c) 2014-2018 Sergio Lucia, Alexandru Tatulea-Codrean
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

# Original do-mpc imports
from casadi import *
import numpy as NP
import core_do_mpc

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

#    The GP Should be trained BEFORE running the MPC (in the respective file).
#    DO NOT train the GP inside the model function, unless you want to re-train 
#    it every time the model runs.


"""
--------------------------------------------------------------------------------
 DEFINE A FUNCTION FOR THE PLASMA MODEL
--------------------------------------------------------------------------------
"""
def model():
    
    """
    --------------------------------------------------------------------------
    template_model: define the non-uncertain parameters
    --------------------------------------------------------------------------
    """
    
    #   We assume that we only have access to a linear model identified
    #   on the glass substrate.
    
    # Glass model (model to which we have access) linearized around Tss = 38 ˚C; Iss = 14 a.u.; Qss = 3 slm; Pss = 3 W
    modelMatG = sio.loadmat('/users/adbonzanini/box sync/berkeley/research/Explicit MPC Paper/Supporting-Data-Files/MIMOmodelGlass')
    A = modelMatG['A']
    B = modelMatG['B']
    C = modelMatG['C']
    



    """
    ----------------------------------------------------------------------------
    template_model: define uncertain parameters, states and controls as symbols
    ----------------------------------------------------------------------------
    """
    # Define the uncertainties as CasADi symbols
    
    p1 = SX.sym("p1")
    p2 = SX.sym("p2")

    # Define the algebraic states as CasADi symbols
    
    T = SX.sym("T")    # Temperature (deviation) [˚C or K]
    I = SX.sym("I")    # Intensity (deviation) [a.u.]
      
    # Define the control inputs as CasADi symbols
    
    q = SX.sym("Q")    # He flowrate [slm]
    P = SX.sym("P")    # Applied power [W]

    # Define time-varying parameters that can change at each step of the prediction and at each sampling time of the MPC controller. For example, future weather predictions

    tv_param_1 = SX.sym("tv_param_1")
    tv_param_2 = SX.sym("tv_param_2")
    
    # Concatenate differential states, algebraic states, control inputs and right-hand-sides
    _x = vertcat(T, I)

    _u = vertcat(q, P)
        
    _p = vertcat(p1, p2)
    
    _z = []

    _tv_p = vertcat(tv_param_1, tv_param_2)

    """
    --------------------------------------------------------------------------
    template_model: define algebraic and differential equations
    --------------------------------------------------------------------------
    """
    '''
    import config
    L = config.L
    if len(config.xopt)>L:
        # Keep the history of the L most recent points
        xopt = config.xopt
        uopt = config.uopt
        
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
        config.w1 = w1
        config.w2 = w2
        
    else:
        # Set GP corrections to zero if no history is available 
        w1 = 0.01
        w2 = 0.01
        config.w1 = w1
        config.w2 = w2
        config.j = config.j+1
    '''
    import config
    w1 = config.w1
    w2 = config.w2
    
    # Concatenate w
    _w = vertcat(w1, w2)

    
    # Define the differential equations
    _xplus = mtimes(A, _x)+mtimes(B, _u) + NP.multiply(_p, _w)
    
    """
    --------------------------------------------------------------------------
    template_model: initial condition and constraints
    --------------------------------------------------------------------------
    """
    # Initial condition for the states
    T0 = 0. # Initial deviation of temperature
    I0 = 0.  # Initial deviation of intensity
    x0 = NP.array([T0, I0])

    # Bounds on the states. Use "inf" for unconstrained states
    T_lb = -5.;			        T_ub = 5.
    I_lb = -15.;			I_ub = 20.
    x_lb = NP.array([T_lb, I_lb])
    x_ub = NP.array([T_ub, I_ub])

    # Bounds on the control inputs. Use "inf" for unconstrained inputs
    q_lb = -1.5;                 q_ub = 7.;
    P_lb = -2.5;                 P_ub = 2.;
    u_lb = NP.array([q_lb, P_lb])
    u_ub = NP.array([q_ub, P_ub])
    u0 = NP.array([3.,3.])

    # Scaling factors for the states and control inputs. Important if the system is ill-conditioned
    x_scaling = NP.array([1.0, 1.0])
    u_scaling = NP.array([1.0, 1.0])

    # Other possibly nonlinear constraints in the form cons(x,u,p) <= cons_ub
    # Define the expresion of the constraint (leave it empty if not necessary)
    cons = vertcat([])
    # Define the lower and upper bounds of the constraint (leave it empty if not necessary)
    cons_ub = NP.array([])

    # Activate if the nonlinear constraints should be implemented as soft constraints
    soft_constraint = 0
    # Penalty term to add in the cost function for the constraints (it should be the same size as cons)
    penalty_term_cons = NP.array([])
    # Maximum violation for the constraints
    maximum_violation = NP.array([0])
    
    '''
    # Define the terminal constraint (leave it empty if not necessary)
    LoadConstraints = sio.loadmat('/users/adbonzanini/box sync/berkeley/research/gaussian process/MIMO GP state feedback substrates/supporting data files/robustConstraints')
    Xf = LoadConstraints['Xtight'][-1][0]
    cons_terminal = Xf[0:2,0:2]
    # Define the lower and upper bounds of the constraint (leave it empty if not necessary)
    cons_terminal_lb = NP.array([-inf]*len(Xf[0:2,:]))
    cons_terminal_ub = NP.array(Xf[0:2,2])
    '''
    
    # Define the terminal constraint (leave it empty if not necessary)
    cons_terminal = vertcat()
    # Define the lower and upper bounds of the constraint (leave it empty if not necessary)
    cons_terminal_lb = NP.array([])
    cons_terminal_ub = NP.array([])


    """
    --------------------------------------------------------------------------
    template_model: steady-states
    --------------------------------------------------------------------------
    """
    
    # Offset-free tracking parameters
    Bd = np.eye(2,2)
    Cd = np.zeros((2,2))
    Haug = np.eye(2,2)
    Aaug = np.vstack([np.hstack([np.eye(2,2)-A, -B]),np.hstack([np.dot(Haug,C), np.zeros((2,2))])])
    
    # Function to calculate the respective inputs
    def spCalculator(yspIn, spYc, Aaug, Bd, Cd): 
    
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
        lbVar +=  [-inf, -inf]
        ubVar +=[inf, inf]
        Var0 +=[0, 0]
        
        Usp = MX.sym('Usp', 2)
        spVar +=[Usp]    
        lbVar +=[-inf, -inf]
        ubVar +=[inf, inf]
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
        feasMeasure = solver.stats()['return_status']
        print('----------------------------------------------')
        print(feasMeasure)
        print('----------------------------------------------')
        
        # Baug = np.vstack([np.dot(Bd, spYc), yspIn-np.dot(np.dot(Haug,Cd), spYc)])
        # spCon = np.linalg.solve(Aaug, Baug)
        # ySP = np.dot(C, spCon[0:2])
        
        # spCon[0:2] = ySP
        # xxsp = np.linalg.solve(C, yspIn-spYc)
        # uu = np.linalg.solve(B, xxsp - np.dot(A, xxsp))
        # spCon = np.array([[float(yspIn[0])], [float(yspIn[1])], [float(uu[0])], [float(uu[1])]])
        
        return sp_opt
    
    # Steady states
    xsp = np.array([[2.], [0.]]).reshape(-1,1)  
    spYc = np.array([[0], [0]]) # No correction due to GP (yet)
    ss = spCalculator(xsp, spYc, Aaug, Bd, Cd)
    
    """
    --------------------------------------------------------------------------
    template_model: cost function
    --------------------------------------------------------------------------
    """
    # Weight matrices
    Q = 1*np.array([[10, 0], [0, 0.1]])
    R = 1*np.array([[5, 0], [0, 2]])
    
    # Define the cost function
    # Lagrange term
    lterm = mtimes(mtimes((_x-ss[0:2,0]).T,Q),(_x-ss[0:2,0]))
    # Mayer term
    mterm =  mtimes(mtimes((_x-ss[0:2,0]).T,Q),(_x-ss[0:2,0])) + mtimes(mtimes((_u-ss[2:,0]).T,R),(_u-ss[2:,0]))
    # Penalty term for the control movements
    rterm = mtimes(mtimes((_u-ss[2:,0]).T,R),(_u-ss[2:,0])) 



    """
    --------------------------------------------------------------------------
    template_model: pass information (not necessary to edit)
    --------------------------------------------------------------------------
    """
    model_dict = {'x':_x,'u': _u, 'rhs':_xplus,'p': _p, 'z':_z,'x0': x0,'x_lb': x_lb,'x_ub': x_ub, 'u0':u0, 'u_lb':u_lb, 'u_ub':u_ub, 'x_scaling':x_scaling, 'u_scaling':u_scaling, 'cons':cons,
    "cons_ub": cons_ub, 'cons_terminal':cons_terminal, 'cons_terminal_lb': cons_terminal_lb,'tv_p':_tv_p, 'cons_terminal_ub':cons_terminal_ub, 'soft_constraint': soft_constraint, 
    'penalty_term_cons': penalty_term_cons, 'maximum_violation': maximum_violation, 'mterm': mterm,'lterm':lterm, 'rterm':rterm}

    model = core_do_mpc.model(model_dict)

    return model

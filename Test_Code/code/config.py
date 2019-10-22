
#   This script is used to store and update global variables that
#   are used between scripts and functions (such as the GP corrections)


#   Move the GP regression in this file and run it inside the templeate_model.py script
#    Try introducing w as a time-varying parameter (tv_param)
import numpy as NP
global w1, w2, j, L
w1 = 99
w2 = 99
j=0
L=1


# xopt = NP.zeros((1,2))
# uopt = NP.zeros((1,2))

# def init():
#     global w1, w2, j, L
#     w1=99
#     w2=99
#     j = 0
#     L = 1
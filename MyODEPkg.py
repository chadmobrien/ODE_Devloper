# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 21:51:56 2014

@author: Chad
"""
from math import ceil
import numpy as np
AbsTol = 1e-3
RelTol = 1e-3


def ODE(diffFunc, tspan, initcond, **Opt):
#    - > diffFunc - the function which defines the first order ODE
#    - > tspan - time Domain
#    - > initcond - initial conditions satisfied at tspan[0]
#    - > **Opt  - Optional EVENT functions and SOLVER 
# Event Function
# detects 0 value event
# value, direction, isterminal = EventFunction(t,y)
#   DIRECTION values: -1 decreasing, 0 both directions, 1 increasing
#   TERMINAL values: 1 yes/stop, 0 for no
#  example for falling object stop when hits ground
#   EVENT  = [y[0]]
#   example when object is at height Y
#  EVENT = [y[0]-Y]
#   EVENT must return Tuple (EVENT, DIRECTION, ISTERMINAL)
    # Flag- stop 

    n = np.size(tspan)        # find the number of iterations
    n1 = np.size(initcond)    # finds nummber of differential equations

    #Initializing storage array for solutions
    T = np.empty((1,1))
    Y = np.empty((1,n1))
    # Initializing 0

    T[0] = tspan[0]
    Y[0] = initcond
    print "Starting Solution..."
    print "T0:" + str(T[0])
    print "Y0:" + str(Y[0])
    
    
    # Setups for Events
    if 'EVENTS' in Opt:
        checkEvents = True
        YE = np.empty((0,n1))
        TE = np.empty((0,1))
        IE = np.empty((0,1))
        EventFunc = Opt['EVENTS']
    else:
        checkEvents = False
    
    # If Specific Solver is given
    if 'SOLVER' in Opt:
        Solver = Opt['SOLVER']
        print "Using AbsTol = " + str(AbsTol)
        print "Using AbsTol = " + str(RelTol)
    else:
        Solver = EulerSolve  # Set Default to Euler Solve in no Input Given

    print "Using Solver: " + str(Solver).split()[1]
    for index in range(1,n):
        T_n = tspan[index]
        # Initial Conditions
        Y_nm1 = Y[index-1]
        T_nm1 = T[index-1]
       
        # Step Setup
        deltaT = T_n - T_nm1
        deltaT = ceil(deltaT*10000)/10000
        

        # Current Step
        Y_n = Solver(diffFunc,T_nm1, Y_nm1, deltaT)
        
        # Appending Solution Vectors
        Y = np.vstack((Y,Y_n))
        T = np.vstack((T,T_n))
        # Check for Events
        if checkEvents:
           iterEvent, iterFlag  = CheckEvents(EventFunc, T,Y, index)
           # checks if Event occurs 
           if iterEvent:
               YE = np.vstack((YE,Y[index]))
               TE = np.vstack((TE,T[index]))
               IE = np.vstack((IE,index))
               # if iterFlag == False, stop sim
               if iterFlag:
                   break;
    if checkEvents:
        outRes  = T,Y,TE,YE,IE
    else:
        outRes = T,Y
        
    return outRes

def EulerSolve(SysDiff,t0, y0, dt):
    dy = SysDiff(t0,y0)
    ynew = y0 + dy*dt
    return ynew
        
def ModifiedEulerSolve(diffFunc,t, initcond, deltaT):
    global AbsTol               # Controls Absolute Error Break
    global RelTol               # Control Relative Error Break
    FullStep = initcond + diffFunc(t,initcond)*deltaT # First Step Approxiation
    ym = initcond + 0.5*deltaT*diffFunc(t,initcond)   # mean value calculation
    tm = t+0.5*deltaT                                 # midpoint time
    y2 = initcond + diffFunc(tm,ym)*deltaT            # second step Approximation
    abs_error = np.sum(y2*y2)-np.sum(FullStep*FullStep)# error calculation
    rel_error = abs_error/np.sum(FullStep*FullStep)     # error calculation
    FullStep = y2                                       # setting looping previous value to compare
    while rel_error > RelTol or abs_error > RelTol:     # loop
        ym = initcond + 0.5*deltaT*diffFunc(t,initcond)
        y2 = initcond + diffFunc(tm,ym)*deltaT
        abs_error = np.sqrt(np.sum(y2*y2)-np.sum(FullStep*FullStep))
        rel_error = abs_error/np.sqrt(np.sum(FullStep*FullStep))
        FullStep = y2
    return y2
     
def CheckEvents(EventFunc, T,Y, c_iter):
    EventOccurs = False
    flag = False
    # YM1 is the previous time step solution
    tnm1 = T[c_iter-1]
    ynm1 = Y[c_iter-1]
    # Current Step
    tn = T[c_iter]
    yn = Y[c_iter]
    value_nm1 = EventFunc(tnm1,ynm1)[0]
    value_n, direction, isterminal = EventFunc(tn,yn)
    #value  = thisEventOut[0]
    #direction = thisEventOut[1]
    #isterminal = thisEventOut[2]
    ubound = np.size(value_n)
    for index in range(0,ubound):
        # Check if Zero is Detected
        if not (value_n[index] == value_nm1[index]):
            if np.allclose(abs(value_n[index]),0.0):
                flag = isterminal[index] # checks if event stops sim
                EventOccurs = True
                break
                # check direction
            elif np.sign(value_n[index]) != np.sign(value_nm1[index]):
                EventOccurs = True
                if direction[index] == 0 : # Catch both decreasing and increasing
                    flag = isterminal[index] #checks if event stops sims
                    break
                    
                elif np.sign((value_n[index]-value_nm1[index])) == direction[index]:
                    flag = isterminal[index] # checks if condition stops sims
                    break
        
    return EventOccurs, flag
    
    
                
                
        
            
        
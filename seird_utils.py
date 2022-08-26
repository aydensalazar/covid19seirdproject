# Utility functions for SEIRD model
from scipy.integrate import odeint
import numpy as np

# global variables
D = 4.0 # CONTAGIOUS PERIOD; DEFAULT: infections lasts four days
gamma = 1.0 / D # Default value

def deriv(y, t, N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end):
    """Calculates the end result for ordinary differential equations for SEIRD values.
    
    Parameters
    ----------
    y : tuple
        Contains current S, E, I, R, D values

    t : int
        Time (days)

    N : int
        Population

    beta : float
        Expected amount of people an infected person infects each day

    gamma : float
        The proportion of infected recovering each day (γ = 1/D)

    delta : float
        Length of incubation period

    alpha : float
        Fatality rate

    rho : float 
        Rate at which individuals die (equals 1/days from infected until death)
    
    R_0_start : float
        Value for R0 on first day

    k : float
        Factor that lets us vary how quickly R_0 declines

    x0 : float
        x-value of the inflection point 

    R_0_end : float
        Value for R0 on last day

    Returns
    ----------
    dSdt : int
        Derivative of S with respect to time
    
    dEdt : int
        Derivative of E with respect to time
    
    dIdt : int 
        Derivative of I with respect to time
        
    dRdt : int
        Derivative of R with respect to time
        
    dDdt : int
        Derivative of D with respect to time
    """
    S, E, I, R, D = y
    dSdt = -beta(t, R_0_start, k, x0, R_0_end) * S * I / N
    dEdt = beta(t, R_0_start, k, x0, R_0_end) * S * I / N - delta * E
    dIdt = delta * E - (1 - alpha) * gamma * I - alpha * rho * I
    dRdt = (1 - alpha) * gamma * I
    dDdt = alpha * rho * I
    return dSdt, dEdt, dIdt, dRdt, dDdt

def logistic_R_0(t, R_0_start, k, x0, R_0_end):
    """Apply a logistic function to R_0 to mimick a continuous change rather than a sharp change in value.
    
    Parameters
    ----------
    t : int
        Time (days)
    
    R_0_start : float
        Value for R0 on first day

    k : float
        Factor that lets us vary how quickly R_0 declines

    x0 : float
        x-value of the inflection point 

    R_0_end : float
        Value for R0 on last day

    Returns 
    --------
    Logistic value of R0 at time step t
    """
    return (R_0_start-R_0_end) / (1 + np.exp(-k*(-t+x0))) + R_0_end

def beta(t, R_0_start, k, x0, R_0_end):
    """The expected amount of people an infected person infects each day
    
    Parameters
    ----------
    t : int
        Time (days)
    
    R_0_start : float
        Value for R0 on first day

    k : float
        Factor that lets us vary how quickly R_0 declines

    x0 : float
        x-value of the inflection point 

    R_0_end : float
        Value for R0 on last day

    Returns 
    --------
    Logistic value of R0 at time step t, with gamma applied. 
    
    """
    return logistic_R_0(t, R_0_start, k, x0, R_0_end) * gamma

def less_than_zero(value, decrement):
    """ Decrements R0 by a specified number. If difference is negative, return 0

    Parameters
    ----------
    value : float
        R0 value 

    decrement : float
        Amount of decay for R0 value

    Returns
    ----------
    int
        New value for R0
    
    """
    if value - decrement < 0:
        return 0
    else:
        return value-decrement
# import necessary libraries
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from seird_utils import *

# configure the application dashboard settings
root = Tk()
root.title('Create your own SEIRD Model for COVID-19')
root.geometry('2000x2000')

def plotseird(t, S, E, I, R, D=None, L=None, R0=None, Alpha=None, CFR=None):
    """ Performs necessary plotting for SEIRD model parameters. 

    Parameters
    ----------
    t : int
        Time (days)
    S : int
        Number of susceptible individuals at time t
    E : int
        Number of exposed individuals at time t
    I : int
        Number of infected individuals at time t
    R : int
        Number of recovered individuals at time t
    D : int
        Number of deceased individuals at time t
    R0 : float
        Basic reproduction number 
    Alpha : float
        Fatality rate
    CFR : float
        Case Fatality Rate (CFR) - the total number of deaths as a proportion of reported cases at time t

    Returns
    -------
    None
        Plots SEIRD Graph

    Code sourced from https://towardsdatascience.com/infectious-disease-modelling-part-i-understanding-sir-28d60e29fdfc  
    """
    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
    ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
    if D is not None:
        ax.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Dead')
        ax.plot(t, S+E+I+R+D, 'c--', alpha=0.7, linewidth=2, label='Total')
    else:
        ax.plot(t, S+E+I+R, 'c--', alpha=0.7, linewidth=2, label='Total')

    ax.set_xlabel('Time (days)')

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend(borderpad=2.0)
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    if L is not None:
        plt.title("Lockdown after {} days".format(L))
    plt.show();

    if R0 is not None or CFR is not None:
        f = plt.figure(figsize=(12,4))

    if R0 is not None:
        # sp1
        ax1 = f.add_subplot(121)
        ax1.plot(t, R0, 'b--', alpha=0.7, linewidth=2, label='R_0')

        ax1.set_xlabel('Time (days)')
        ax1.title.set_text('R_0 over time')
        
        ax1.yaxis.set_tick_params(length=0)
        ax1.xaxis.set_tick_params(length=0)
        ax1.grid(b=True, which='major', c='w', lw=2, ls='-')
        legend = ax1.legend()
        legend.get_frame().set_alpha(0.5)
        for spine in ('top', 'right', 'bottom', 'left'):
            ax.spines[spine].set_visible(False)

    if Alpha is not None:
        # sp2
        ax2 = f.add_subplot(122)
        ax2.plot(t, Alpha, 'r--', alpha=0.7, linewidth=2, label='alpha')

        ax2.set_xlabel('Time (days)')
        ax2.title.set_text('fatality rate over time')
        
        ax2.yaxis.set_tick_params(length=0)
        ax2.xaxis.set_tick_params(length=0)
        ax2.grid(b=True, which='major', c='w', lw=2, ls='-')
        legend = ax2.legend()
        legend.get_frame().set_alpha(0.5)
        for spine in ('top', 'right', 'bottom', 'left'):
            ax.spines[spine].set_visible(False)

        plt.show()

# setting global variables 
D = 4.0 # CONTAGIOUS PERIOD; DEFAULT: infections lasts four days
gamma = 1.0 / D # Default value

# adding labels to user interface
Label(root, text='Population Size:', font=(None, 30)).grid(row=0, column=0)
Label(root, text='Mortality Rate:', font=(None, 30)).grid(row=1, column=0)
Label(root, text='Contagious Period (Days):', font=(None, 30)).grid(row=2, column=0)
Label(root, text='Incubation Time (Days):', font=(None, 30)).grid(row=3, column=0)
Label(root, text='Days from Infection to Death:', font=(None, 30)).grid(row=4, column=0)
Label(root, text='Select policy/policies to implement for a COVID-19 Response:', font=(None, 15)).grid(row=5, column=0)

# inserting text boxes for user input

populationinput = Entry(root, width=75)
populationinput.insert(0, 1000000)
mortalityrateinput = Scale(root, from_=0, to=1, orient=HORIZONTAL, resolution=0.1, length=500, width = 30)
mortalityrateinput.set(0.3)
contagiousperiodinput = Entry(root, width=75)
contagiousperiodinput.insert(0, 4)
incubationtimeinput = Entry(root, width=75)
incubationtimeinput.insert(0, 5)
infectiontodeathinput = Entry(root, width=75)
infectiontodeathinput.insert(0, 9)

businesses_var = IntVar()
businesses = Checkbutton(root, text="Close non-essential businesses", variable =businesses_var)
schools_var = IntVar()
schools = Checkbutton(root, text="Close schools and universities", variable =schools_var)
quarantine_var = IntVar()
quarantine = Checkbutton(root, text="Stay at home orders (quarantine)", variable =quarantine_var)
testing_var = IntVar()
testing = Checkbutton(root, text="Aggressive testing", variable =testing_var)
travel_var = IntVar()
travel = Checkbutton(root, text="Travel ban", variable =travel_var)

# configuring arrangement of text boxes for user interface

populationinput.grid(row=0, column=1)
mortalityrateinput.grid(row=1, column=1)
contagiousperiodinput.grid(row=2, column=1)
incubationtimeinput.grid(row=3, column=1)
infectiontodeathinput.grid(row=4, column=1)
businesses.grid(row=6, column=0)
schools.grid(row=7,column=0)
quarantine.grid(row=8,column=0)
testing.grid(row=9, column=0)
travel.grid(row=10, column=0)


def buildgraph():
    """ Builds a SEIRD model based on the inputs given in the tkinter entries above 

    Parameters
    -------
    None

    Returns
    -------
    None
        Calls function to plots SEIRD Graph

    """

    # grab the necessary user inputs from the text boxes
    N = float(populationinput.get())
    alpha = float(mortalityrateinput.get())
    D = float(contagiousperiodinput.get())
    gamma = 1.0 / D
    delta = 1/float(incubationtimeinput.get())
    rho = 1/float(infectiontodeathinput.get())
    S0, E0, I0, R0, D0 = N-1, 1, 0, 0, 0
    t = np.linspace(0, 99, 100) # Grid of time points (in days)
    y0 = S0, E0, I0, R0, D0 # Initial conditions vector
    R_0_start, k, x0, R_0_end = 6.0, 1, 30, 6.0 # initial R0 Values
    policy_implemented = False

    # apply special conditions, if specified by user
    if businesses_var.get():
        policy_implemented = True
        R_0_end = less_than_zero(R_0_end, 0.5)

        # Integrate the SIR equations over the time grid, t:
        ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end))

        S, E, I, R, D = ret.T
        R0_over_time = [logistic_R_0(i, R_0_start, k, x0, R_0_end) for i in range(len(t))]

    if schools_var.get():
        policy_implemented = True
        R_0_end = less_than_zero(R_0_end, 0.1)

        # Integrate the SIR equations over the time grid, t:
        ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end))

        S, E, I, R, D = ret.T
        R0_over_time = [logistic_R_0(i, R_0_start, k, x0, R_0_end) for i in range(len(t))]

    if quarantine_var.get():
        policy_implemented = True
        R_0_end = less_than_zero(R_0_end, 0.4)

        # Integrate the SIR equations over the time grid, t:
        ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end))

        S, E, I, R, D = ret.T
        R0_over_time = [logistic_R_0(i, R_0_start, k, x0, R_0_end) for i in range(len(t))]

    if testing_var.get():
        policy_implemented = True
        R_0_end = less_than_zero(R_0_end, 0.5)

        # Integrate the SIR equations over the time grid, t:
        ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end))

        S, E, I, R, D = ret.T
        R0_over_time = [logistic_R_0(i, R_0_start, k, x0, R_0_end) for i in range(len(t))]

    if travel_var.get():
        policy_implemented = True
        R_0_end = less_than_zero(R_0_end, 0.2)

        # Integrate the SIR equations over the time grid, t:
        ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end))

        S, E, I, R, D = ret.T
        R0_over_time = [logistic_R_0(i, R_0_start, k, x0, R_0_end) for i in range(len(t))]

    if not policy_implemented:
        R_0_start, k, x0, R_0_end = 6.0, 1, 70, 6.0

        # Integrate the SIR equations over the time grid, t:
        ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha, rho, R_0_start, k, x0, R_0_end))

        S, E, I, R, D = ret.T
        R0_over_time = [logistic_R_0(i, R_0_start, k, x0, R_0_end) for i in range(len(t))]

    # call plotting function 
    plotseird(t, S, E, I, R, D, R0=R0_over_time)
    

GenerateGraphButton = Button(root, text='Generate Graph', command = buildgraph)
GenerateGraphButton.grid(row=15, column=1)

root.mainloop()

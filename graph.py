#module for plotting figures for the protein class in Main.py
import numpy as np
import matplotlib.pyplot as plt


def generatechargeplot(sequences):
    chargeplot = plt.figure()
    for i in range(0,int(len(sequences))): #generates the y values
        sequences[i].chargeplot() 

    plt.grid() # adds gridlines

    plt.xlabel('pH') 
    plt.ylabel('Relative Charge')
    plt.legend() 

    #configure axes
    ax = plt.gca() #get current axes (then spines to configure axes)
    ax.spines['top'].set_color('none') #removes top black line on plot
    ax.spines['bottom'].set_position('zero') #ensures x axis always crosses at y = 0
    ax.spines['right'].set_color('none') #removes right black line on plot

    #x axis ticks  
    start,end = ax.get_xlim() #gets min and max of x axis
    end = end+1
    ax.set_xticks(np.arange(0,round(end),1)) #sets major ticks on x axis
    ax.set_xticks(np.arange(round(start),round(end),0.25),True) #sets minor ticks on x axis

    #y axis ticks 
    start,end = ax.get_ylim() #gets min and max of y axis
    end = end 
    ax.set_yticks(np.arange(round(start),round(end),1),True) #sets minor ticks on y axis
    
    return chargeplot

def generatemassplot(sequences):
    massplot = plt.figure()
    for i in range(0,int(len(sequences))): #generates the y values
        sequences[i].massplot()

    plt.xlabel('Protein')
    plt.ylabel('Mass (kDa)')
    plt.legend()

    #configure axes
    bx = plt.gca() #get current axes (then spines to configure axes)
    bx.spines['top'].set_color('none') #removes top black line on plot
    bx.spines['bottom'].set_position('zero') #ensures x axis always crosses at y = 0
    bx.spines['right'].set_color('none') #removes right black line on plot
    bx.set_axisbelow(True)

    #y axis ticks 
    start,end = bx.get_ylim() #gets min and max of y axis
    start = 0
    end = end + 1
    bx.set_yticks(np.arange(0,round(end),5)) #sets major ticks on x axis
    bx.set_yticks(np.arange(round(start),round(end),1),True) #sets minor ticks on y axis
    bx.yaxis.grid(which='major')

    return massplot
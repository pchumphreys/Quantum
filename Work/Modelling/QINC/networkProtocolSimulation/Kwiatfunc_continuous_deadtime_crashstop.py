"""
numerical simulation of the 'MidpointSource' protocol as described in Jones et al (NJP 18 (2016) 083015)
"""

import numpy as np
import random

def Kwiat_crashStop(n_iterations=3000,pfc = 0.3,pout=0.3,pm=0.1,d=50.0, Eg=1.):
    """
    The simulation is done according to fixed rounds between which communication takes place, as in Jones et al. 
    Room for improvement is to do the scheme in a continuous manner
    TODO: do not for every emission event invoke a random event. do it smarter; ths is too slow
    """
    c=0.2

    n=n_iterations #no of times we repeat the experiment

    Ploss = 10.**(-0.2*(d/20.)) #losses over the distance of d/2

    pr = 0.5*Ploss*pfc*pout #probability that photon from entangled pair is latched in right repeater
    pl = 0.5*Ploss*pfc*pout #" " "in left repeater
    p = pl*pm*pr   #total probability of succesful distant entanglement generation -> not used now (only indirectly)
    
    num_random = 1000
    rand_l = np.random.uniform(0,1,num_random) < pl
    rand_r = np.random.uniform(0,1,num_random) < pr
    rand_m = np.random.uniform(0,1,num_random) < pm
    x_l = 0
    x_r = 0
    x_m = 0

    time_steps_to_comm = d/(c*Eg)
    
    data = np.zeros(n);
    
    for i in range(0, n):
        time = 0
        succ = 0
        t_stored_A = 0
        t_stored_B = 0
        a_latched = False
        b_latched = False

        while succ == 0:

            if rand_m.item(x_m): #probability that there was a succesful entangled state sent
                    if not(a_latched):
                        if rand_l.item(x_l): #probability of succesful latching at Alice (pl)
                            t_stored_A = time
                            a_latched = True
                            
                        x_l += 1
                        if x_l == num_random:
                            rand_l = np.random.uniform(0,1,num_random) < pl
                            x_l = 0

                    if not(b_latched):
                        if rand_r.item(x_r):
                            t_stored_B = time
                            b_latched = True

                        x_r += 1
                        if x_r == num_random:
                            rand_r = np.random.uniform(0,1,num_random) < pr
                            x_r = 0

            x_m += 1
            if x_m == num_random:
                rand_m = np.random.uniform(0,1,num_random) < pm
                x_m = 0


            if a_latched and time >= (t_stored_A + time_steps_to_comm): #there is information available about the oldest memory
                
                if b_latched and  t_stored_A == t_stored_B:
                    succ += 1
                else:
                    a_latched = False

            if b_latched and time >= (t_stored_B + time_steps_to_comm): #there is information available about the oldest memory
                
                b_latched = False
                            
            
            time += 1
            
        data[i] = time

        if i%200==0:
            print '%d out of %d done'%(i+1,n)

    data = data*Eg
    output = np.mean(data)
    std_output = np.std(data)

    return [output, std_output]




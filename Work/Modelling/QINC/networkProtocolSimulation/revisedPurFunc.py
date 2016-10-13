# -*- coding: utf-8 -*-
"""
numerical simulation of the 'MidpointSource' protocol as described in Jones et al (NJP 18 (2016) 083015)
"""

import numpy as np
import random

def Pur(NoMem,n_iterations=3000,pfc = 0.3,pout=0.3,d=50.0,fid_bound =0.51, Eg=1.,Sg=20.,Pg=10.,mLt=400.):
    """
    The simulation is done according to fixed rounds between which communication takes place, as in Jones et al. 
    Room for improvement is to do the scheme in a continuous manner
    TODO: do not for every emission event invoke a random event. do it smarter; ths is too slow
    """
    c=0.2

    Ploss = 10.**(-0.2*(d/20.)) #losses over the distance of d/2
 
    P = Ploss * pout * pfc
    
    pSucc = P * (1 - 0.5 * P) # Prob of click 

    max_noDec = int(-mLt*np.log(2*fid_bound -1 )) #if this is exceeded, stop trying to attempt generating entanglement

    data = np.zeros((n_iterations,3))

    Mem_NoDecs = np.zeros(NoMem) # Track decoherence
    t_avail = np.zeros(NoMem) # Track the time when states will be available
    t_purify_known = np.zeros(NoMem) # Track the time when purification results will be known

    # Pregenerate lots of random numbers for speed
    num_rand = 1000
    successes = np.random.uniform(0, 1, num_rand) < pSucc
    x = 0
   
    for i in range(n_iterations):

        # Indices to track where to store stuff, and what is going on
        num_stored = 0
        num_mem_avail = NoMem
        num_ent_created = 0
        purifications_in_prog = 0
        t_avail_zero_ind = 0
        t_purify_known_zero_ind = 0
        Mem_NoDecs_avail_zero_ind = 0

        t = 0 # Time

        #--------------------------------
        # Loop until success
        while True:

            #--------------------------------
            # Check for available events
            # If entangled states have made it to the beamsplitter, and info has propagated back
            if num_ent_created and t >= t_avail.item(t_avail_zero_ind):

                t_avail_zero_ind = (t_avail_zero_ind + 1) % NoMem
                num_ent_created -= 1 # Reduce index by one

                #--------------------------------
                # Check for successful events
                
                succ_this_round = successes.item(x)
                x += 1
                if x == num_rand: # Need to make new random numbers!
                    successes = np.random.uniform(0, 1, num_rand) < pSucc
                    x = 0

                if succ_this_round:  #Info that this entanglement round succeeded

                    Mem_NoDecs_next_ind_to_write = (num_stored + Mem_NoDecs_avail_zero_ind) % NoMem 
                    Mem_NoDecs.itemset(Mem_NoDecs_next_ind_to_write, num_ent_created)  # This gives the number of entanglement gen events that happened after this entangled state was created 
                    num_stored += 1

                    if num_stored == 2: # Two memories are enough to purify!
               
                        t += Pg
               
                        # Check now whether purify will succeed, currently we only consider the possibility of one purification succeeding at a time
                        if np.random.uniform(0, 1) < (1/8.0): 
                            Mem_NoDecs_first_two_inds = (np.arange(0,2) + Mem_NoDecs_avail_zero_ind) % NoMem
                            succ_Decs = np.sum(Mem_NoDecs[Mem_NoDecs_first_two_inds])
                            t_succ = t + d/c
                            break # In this case, can just break. In principle could keep going for more complicated routines

                        else:
                            t_purify_known_next_ind_to_write = (purifications_in_prog + t_purify_known_zero_ind) % NoMem 
                            t_purify_known.itemset(t_purify_known_next_ind_to_write, t + d/c)
                            purifications_in_prog += 1
                            Mem_NoDecs_avail_zero_ind = (Mem_NoDecs_avail_zero_ind + 2) % NoMem
                            num_stored -= 2
                            num_mem_avail += 1
                else:
                    num_mem_avail += 1

            #--------------------------------
            # Check if unsuccessful purification info has made it back
                
            if purifications_in_prog and t >= t_purify_known.item(t_purify_known_zero_ind): # Check whether the purification result is known yet!
                        
                    purifications_in_prog -= 1
                    t_purify_known_zero_ind = (t_purify_known_zero_ind + 1) % NoMem
                
                    num_mem_avail += 1            
            
            #--------------------------------
            # Check if can make entanglement using available memories

            if num_mem_avail and not(purifications_in_prog): # Only create entangled states when purification isnt in progress
                    t_avail_next_ind_to_write = (t_avail_zero_ind + num_ent_created) % NoMem 
                    t_avail.itemset(t_avail_next_ind_to_write, t + Eg + d/c)
                    num_ent_created += 1
                    num_mem_avail -= 1

                    t += Eg + Sg

                    # Induce decoherence from this entanglement generation
                    Mem_NoDecs += 1
                    
                    # Check for dead mem
                    if num_stored and Mem_NoDecs.item(Mem_NoDecs_avail_zero_ind) > max_noDec:
                        Mem_NoDecs_avail_zero_ind = (Mem_NoDecs_avail_zero_ind + 1) % NoMem
                        num_stored -= 1
                        num_mem_avail += 1

            else: # If cannot do anything else, need to fast forward
                if purifications_in_prog and ((num_ent_created and t_purify_known.item(t_purify_known_zero_ind) < t_avail.item(t_avail_zero_ind)) or not(num_ent_created)):
                    t = t_purify_known.item(t_purify_known_zero_ind)
                elif num_ent_created:
                    t = t_avail.item(t_avail_zero_ind)
                else: # This really shouldnt ever happen!
                    print "System seems to be stuck!"
                    # print "time: ", t
                    # print "t_avail: ", t_avail
                    # print "num_ent_created: ", num_ent_created

       
        F = (1. + np.exp( -(succ_Decs)/mLt ))/2 #fidelity decays with this number of attempts after the successful run 

        data[i,:] = [t_succ, F,succ_Decs]

        if i%200==0:
            print '%d out of %d done'%(i+1,n_iterations)

    output = np.mean(data, axis=0)
    std_output = np.std(data,axis=0)

    return np.append(output,std_output) #output: [duration, fidelity]




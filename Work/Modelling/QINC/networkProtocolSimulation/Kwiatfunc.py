"""
numerical simulation of the 'MidpointSource' protocol as described in Jones et al (NJP 18 (2016) 083015)
"""

import numpy as np
import random

def Kwiat(NoMem,n_iterations=3000,K=None,pfc = 0.3,pout=0.3,pm=0.1,d=50.0):
    """
    The simulation is done according to fixed rounds between which communication takes place, as in Jones et al. 
    Room for improvement is to do the scheme in a continuous manner
    TODO: do not for every emission event invoke a random event. do it smarter; ths is too slow
    """

    Eg = 1. 
    Sg = 20. 
    mLt = 400.
    c=0.2

    n=n_iterations #no of times we repeat the experiment

    Ploss = 10.**(-0.2*(d/20.)) #losses over the distance of d/2

    pr = 0.5*Ploss*pfc*pout #probability that photon from entangled pair is latched in right repeater
    pl = 0.5*Ploss*pfc*pout #" " "in left repeater
    p = pl*pm*pr   #total probability of succesful distant entanglement generation -> not used now (only indirectly)
    
    #print 'losses over distance %d km are %.3f'%(d,Ploss)
    #print 'pm',pm
    #print 'pr = pl =',pr
    #print 'K: %d'%K

    if K==None:
        K = int(NoMem/(pl*pr))# number of latch attempts per round. set to optimum according to Jones et al. 
                #-> this will not for us the optimum since we have decoherence induced after succesful latching. 
                # we could already stop trying to entangle after all the memories are full though.

    round_time = K*(Eg+Sg)+d/c

    data = np.zeros((n,2));
    for i in range(0, n):
        time = 0
        succ = 0
        #memA = 0
        #memB = 0
        rounds = 0
        while succ == 0:
            j_A = np.array([]) #the bins in which a photon is latched at A
            j_B = np.array([]) # " " at B
            time = time + round_time #since we work in discrete rounds, the time always increases by this fixed number.

            for j in np.arange(K):#there are K emission events at the source
                if random.uniform(0, 1) < pm: #probability that there was a succesful entangled state sent
                    if (len(j_A)<NoMem):
                        if random.uniform(0,1) < pr: #probability of succesful latching at Alice (pr)
                            j_A = np.append(j_A,j) #fill a memory of Alice
                            #print 'appending to ja'
                    if (len(j_B)<NoMem):
                        if random.uniform(0,1) < pl:
                            j_B = np.append(j_B,j) #fill a memory of Bob
                            #print 'appending to jb'
                    if (len(j_A)==NoMem) or (len(j_B)==NoMem) and succ == 0:#one of the memories is already full!
                        j_success = np.intersect1d(j_A,j_B) #only taken into account the actual size of the memory
                        #print 'j_a',j_A
                        #print 'j_b',j_B
                        #print j_success
                        succ += len(j_success) #we increment success by the number of succesful entanglement events.
                        #print succ
                    if (len(j_A)==NoMem) and (len(j_B)==NoMem):#both memories are already full!
                        break #we can break out of the for loop for this round now, since both the memories are full.
                        # we have to give both memories the chance to fill completely before breaking out of the loop to know 
                        # how much decoherence due to entanglement generation we have on both sidews.
                       
            rounds += 1        

        if succ >1: #this probability is very low. right now do nothing with it.
            pass

        #print 'rounds needed:%d'%rounds
        #print j_success
        #print len(j_A)
        #print len(j_B)
        #determine when Alice resp Bob stopped attempting to generate entanglement (memory full, or round stopped)
        if len(j_A)<NoMem:
            j_stopA = K
        else:
            j_stopA = j_A[-1]
        if len(j_B)<NoMem:
            j_stopB = K
        else:
            j_stopB = j_B[-1]
        #number of entanglement events after the succesful bin. this decoherece the entangled state
        NoDecA = j_stopA - j_success[0] 
        NoDecB = j_stopB - j_success[0]
        #print NoDecA, NoDecB

        F = (1. + np.exp( -(NoDecA+NoDecB )/(2*mLt) ))/2 #fidelity decays with this number of attempts after the successful run 
            

        data[i,:] = [time, F]

        if i%200==0:
            print '%d out of %d done'%(i+1,n)

    output = np.mean(data, axis=0)

    return [output[0], output[1]] #output: [duration, fidelity]




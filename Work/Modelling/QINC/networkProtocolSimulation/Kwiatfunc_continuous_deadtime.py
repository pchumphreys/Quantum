"""
numerical simulation of the 'MidpointSource' protocol as described in Jones et al (NJP 18 (2016) 083015)
"""

import numpy as np
import random

def Kwiat(NoMem,n_iterations=3000,pfc = 0.3,pout=0.3,pm=0.1,d=50.0,fid_bound =0.51, Eg=1.,Sg=20.,mLt=400.):
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
    
    #print 'losses over distance %d km are %.3f'%(d,Ploss)
    #print 'pm',pm
    #print 'pr = pl =',pr
    #print 'K: %d'%K

    #fid_bound = (1. + np.exp( -(Max_no_Dec)/(2*mLt) ))/2
    #print fid_bound
    #print np.log(fid_bound - 1/2.)
    max_noDec = int(-mLt*np.log(2*fid_bound -1 )) #if this is exceeded, stop trying to attempt generating entanglement
    print 'max number of attempts after success event:',max_noDec

    #print (1. + np.exp( -(max_noDec)/(mLt) ))/2

    data = np.zeros((n,4));
    for i in range(0, n):
        NoDecs_A = np.zeros(NoMem)
        NoDecs_B = np.zeros(NoMem)
        j_A = np.array([]) #the bins in which a photon is latched at A
        j_B = np.array([]) # " " at B
        dead_time_A = False
        dead_time_B = False
        stop_A = False
        stop_B = False
        time = 0
        succ = 0
        j=0
        #memA = 0
        #memB = 0
        rounds = 0
        while succ == 0:

            time = time + (Eg)  #since we work in discrete rounds, the time always increases by this fixed number.

            if len(j_A)>0 and ((j_A[0] + int(Sg/Eg) + max_noDec) < j):#check if the fidelity bound is reached
                #print j,'stopping entanglement attempts A'
                stop_A = True

            if not (dead_time_A  or stop_A):
                NoDecs_A[:len(j_A)]+=1 #decoherence doesn't depend on successful sending of entangled state
                #print NoDecs_A[:len(j_A)]

            if len(j_B)>0 and ((j_B[0] + int(Sg/Eg)+ max_noDec) < j):#check if the fidelity bound is reached
                #print j,'stopping entanglement attempts B'
                stop_B = True

            if not (dead_time_B or stop_B): #only if running, have a decoherence event.
                #print NoDecs_B[:len(j_B)]
                NoDecs_B[:len(j_B)]+=1

            if random.uniform(0, 1) < pm: #probability that there was a succesful entangled state sent
                if (len(j_A)<NoMem) and dead_time_A == False and stop_A == False:
                    if random.uniform(0,1) < pl: #probability of succesful latching at Alice (pl)
                        j_A = np.append(j_A,j) #fill a memory of Alice
                        dead_time_A = True
                        #print j, 'A dead',j_A
                                                #print 'appending to ja',j
                if (len(j_B)<NoMem) and dead_time_B == False and stop_B == False:
                    if random.uniform(0,1) < pr:
                        j_B = np.append(j_B,j) #fill a memory of Bob
                        dead_time_B = True
                        #print j, 'B dead'
                        #print 'appending to jb',j

            if dead_time_A:
                if (j_A[-1] + Sg/Eg < j): #has the swapping finished yet?
                    #print j,'A back alive'
                    dead_time_A = False

            if dead_time_B:
                if (j_B[-1] + Sg/Eg < j):
                    #print j,'B back alive'
                    dead_time_B = False


            #print j_A
            #print j_B
            # if len(j_A)>0:
            #     print (j_A[0] + d/c/(Eg+Sg))
            #     print j

            if (not dead_time_A) and (not dead_time_B):
                #print 'both alive'
                if (len(j_A)>0): #we only have to check the memory if there is anything saved
                    if ((j_A[0] + d/c/(Eg)) < j): #there is information available about the oldest memory
                        #print 'info available about', int(j_A[0])
                        stop_A = False #we can start generating entanglement again
                        if j_A[0] in j_B: #success!
                            j_success = j_A[0]
                            #print 'success',j_success
                            succ += 1
                        else: #reuse this memory
                            #print 'curretn ja:',j_A
                            #print 'curretn nodecs A',NoDecs_A
                            j_A = np.delete(j_A, 0)
                            NoDecs_A = np.delete(NoDecs_A,0)#also remove the decoherence info for this memory
                            NoDecs_A = np.append(NoDecs_A,0)#add room for decoherence info about new memory to NoDecs
                            #print 'deleted:',j_A
                            #print 'NoDecs new:',NoDecs_A
                if len(j_B)>0:
                    if (j_B[0] + d/c/(Eg)) < j : #we have info about succesful latching in B for j_A
                        stop_B = False
                        #print 'info available about', int(j_B[0])
                        if succ ==0: #we already know whether both a and b latched a photon from the same bin (if statement above)
                            j_B = np.delete(j_B, 0)
                            NoDecs_B = np.delete(NoDecs_B,0)#also remove the decoherence info for this memory
                            NoDecs_B = np.append(NoDecs_B,0)#add room for decoherence info about new memory to NoDecs

            j+=1
            #if j>10000:#use break statement for testing
            #   break   

        if succ >1: #this probability is very low. right now do nothing with it.
            pass

        #print 'success mem A,B=',success_ix_A[0],success_ix_B[0] #check:these should always be indices 0!


        NoDecA = NoDecs_A[0] 
        NoDecB = NoDecs_B[0]
        #print NoDecA, NoDecB

        F = (1. + np.exp( -(NoDecA +NoDecB  )/(2*mLt) ))/2 #fidelity decays with this number of attempts after the successful run 
            

        data[i,:] = [time, F,NoDecA,NoDecB]

        if i%200==0:
            print '%d out of %d done'%(i+1,n)

    output = np.mean(data, axis=0)
    std_output = np.std(data,axis=0)

    return [output[0], output[1],output[2],output[3],std_output[0] ,std_output[1],std_output[2],std_output[3]] #output: [duration, fidelity]




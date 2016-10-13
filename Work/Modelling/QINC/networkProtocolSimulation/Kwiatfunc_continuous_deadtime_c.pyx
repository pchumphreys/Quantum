# cython: profile=True
"""
numerical simulation of the 'MidpointSource' protocol as described in Jones et al (NJP 18 (2016) 083015)
"""

import numpy as np
cimport numpy as np
from libc.stdlib cimport rand, RAND_MAX
from cpython cimport bool
import cython
@cython.boundscheck(False)

def Kwiat(int NoMem,int n_iterations=3000,double pfc = 0.3,double pout=0.3,double pm=0.1,double d=50.0,double fid_bound =0.51, double Eg=1.,double Sg=20.,double mLt=400.):
    """
    The simulation is done according to fixed rounds between which communication takes place, as in Jones et al. 
    Room for improvement is to do the scheme in a continuous manner
    TODO: do not for every emission event invoke a random event. do it smarter; ths is too slow
    """
    cdef double c, Ploss, pr, pl, p, time, F
    cdef int max_noDec, mems_used_A, mems_used_B, succ, j, rounds, ratioEgSg, jA0, jB0
    cdef bool dead_time_A, dead_time_B, stop_A, stop_B

    c=0.2

    Ploss = 10.**(-0.2*(d/20.)) #losses over the distance of d/2

    pr = 0.5*Ploss*pfc*pout #probability that photon from entangled pair is latched in right repeater
    pl = 0.5*Ploss*pfc*pout #" " "in left repeater
    p = pl*pm*pr   #total probability of succesful distant entanglement generation -> not used now (only indirectly)
    
    cdef np.ndarray[np.int_t,
                    ndim=1,
                    negative_indices=False,
                    mode='c'] j_A = np.zeros(NoMem)

    cdef np.ndarray[np.int_t,
                    ndim=1,
                    negative_indices=False,
                    mode='c'] j_B = np.zeros(NoMem)

    cdef np.ndarray[np.int_t,
                    ndim=1,
                    negative_indices=False,
                    mode='c'] NoDecs_A = np.zeros(NoMem)

    cdef np.ndarray[np.int_t,
                    ndim=1,
                    negative_indices=False,
                    mode='c'] NoDecs_B = np.zeros(NoMem)

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

    ratioEgSg = int(Sg/Eg)

    data = np.zeros((n_iterations,4));
    for i in range(0, n_iterations):
        NoDecs_A = np.zeros(NoMem)
        NoDecs_B = np.zeros(NoMem)
        j_A = np.zeros(NoMem) #the bins in which a photon is latched at A
        j_B = np.zeros(NoMem) # " " at B
        mems_used_A = 0
        mems_used_B = 0
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

            if mems_used_A > 0 and ((jA0 + ratioEgSg + max_noDec) < j):#check if the fidelity bound is reached
                #print j,'stopping entanglement attempts A'
                stop_A = True

            if not (dead_time_A  or stop_A):
                NoDecs_A[:mems_used_A]+=1 #decoherence doesn't depend on successful sending of entangled state
                #print NoDecs_A[:len(j_A)]

            if mems_used_B>0 and ((jB0 + ratioEgSg+ max_noDec) < j):#check if the fidelity bound is reached
                #print j,'stopping entanglement attempts B'
                stop_B = True

            if not (dead_time_B or stop_B): #only if running, have a decoherence event.
                #print NoDecs_B[:len(j_B)]
                NoDecs_B[:mems_used_B]+=1

            if rand()/(RAND_MAX*1.0) < pm: #probability that there was a succesful entangled state sent
                if (mems_used_A<NoMem) and dead_time_A == False and stop_A == False:
                    if rand()/(RAND_MAX*1.0) < pl: #probability of succesful latching at Alice (pl)
                        j_A[mems_used_A] = j #fill a memory of Alice
                        jA0 = j_A[0]
                        mems_used_A += 1
                        dead_time_A = True
                        #print j, 'A dead',j_A
                                                #print 'appending to ja',j
                if (mems_used_B<NoMem) and dead_time_B == False and stop_B == False:
                    if rand()/(RAND_MAX*1.0) < pr:
                        j_B[mems_used_B] = j #fill a memory of Alice
                        jB0 = j_B[0]
                        mems_used_B += 1
                        dead_time_B = True
                        #print j, 'B dead'
                        #print 'appending to jb',j

            if dead_time_A:
                if (j_A[mems_used_A-1] + ratioEgSg < j): #has the swapping finished yet?
                    #print j,'A back alive'
                    dead_time_A = False

            if dead_time_B:
                if (j_B[mems_used_B-1] + ratioEgSg < j):
                    #print j,'B back alive'
                    dead_time_B = False


            #print j_A
            #print j_B
            # if len(j_A)>0:
            #     print (j_A[0] + d/c/(Eg+Sg))
            #     print j

            if (not dead_time_A) and (not dead_time_B):
                #print 'both alive'
                if (mems_used_A>0): #we only have to check the memory if there is anything saved
                    if ((jA0 + d/c/(Eg)) < j): #there is information available about the oldest memory
                        #print 'info available about', int(j_A[0])
                        stop_A = False #we can start generating entanglement again
                        if jA0 in j_B: #success!
                            j_success = jA0
                            #print 'success',j_success
                            succ += 1
                        else: #reuse this memory
                            #print 'curretn ja:',j_A
                            #print 'curretn nodecs A',NoDecs_A
                            j_A = np.roll(j_A, -1)
                            j_A[-1] = 0
                            jA0 = j_A[0]
                            NoDecs_A = np.roll(NoDecs_A, -1)
                            NoDecs_A[-1] = 0
                            mems_used_A -= 1
                if mems_used_B>0:
                    if (jB0 + d/c/(Eg)) < j : #we have info about succesful latching in B for j_A
                        stop_B = False
                        #print 'info available about', int(j_B[0])
                        if succ ==0: #we already know whether both a and b latched a photon from the same bin (if statement above)
                            j_B = np.roll(j_B, -1)
                            j_B[-1] = 0
                            jB0 = j_B[0]
                            NoDecs_B = np.roll(NoDecs_B, -1)
                            NoDecs_B[-1] = 0
                            mems_used_B -= 1

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
            print '%d out of %d done'%(i+1,n_iterations)

    output = np.mean(data, axis=0)
    std_output = np.std(data,axis=0)

    return [output[0], output[1],output[2],output[3],std_output[0] ,std_output[1],std_output[2],std_output[3]] #output: [duration, fidelity]




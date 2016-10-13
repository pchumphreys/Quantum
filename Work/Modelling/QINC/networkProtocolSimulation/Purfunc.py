import random
import numpy
#from matplotlib import pyplot as plt
# Simulation of the multiplexed Purification scheme.
# Here the number of carbon memories and the cutoff-fidelity (x) are the inputs.

def Pur(NoMem, n_iterations = 3000, pfc = 0.3,pout=0.3,d=50.0,x = 0, Eg=1.,Sg=20.,mLt=400.):
	pass
	pfc=0.3
	pout = 0.3
	Eg = 1.
	Sg = 20.
	mLt = 400.
	c=0.2
	timepur = 10


	# pfc=0.1
	# pout = 0.1
	# Eg = 2.
	# Sg = 200.
	# mLt = 200.
	# c=0.2
	# d=50.0
	# NoMem = 5
	# timepur = 100

	if (NoMem-2) * (Eg + Sg) >= d/c: #If we have too many memories such that the heralding signal comes back before we would use all memories, then we will use only those first memories that can be all used before the signal comes back
		NoMem = numpy.ceil( d/(c*(Eg+Sg)) + 1 ) #Note that here we allow one more memory that in BK because after generating the first copy we can store it away and use remaining memories for getting a new click
		print " Too many memories. We will only use ", NoMem,  " memories "

	
	Ploss = 10.**(-0.2*(d/20.)) #losses over the distance of d/2

	P = Ploss * pout * pfc
	data = numpy.zeros((n_iterations,2));
	
	for i in range(0, n_iterations):
		timefinal = 0
		succpur=0
		succ=0
		mem = 0
		time3succ=0

		while succpur == 0:



			fid = 1 #we start from fidelity of 1
			time = 0 #time until first click
			rounds = 1 #a round lasts until we try once with all the memories
						
			while succ == 0:
				time = time + Eg + Sg  #progress a trial in time
				mem += 1.
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					succ += 1.
					click = 1. # this will be relevant later, if click=0, this means that we have some entangled state copy left after an unsuccesfull purification attempt
					timesucc = time # time at which we succeeded
				if mem == min(NoMem, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded
					time = max( rounds * (d/c), rounds * mem * (Eg + Sg) ) #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c or a little bit longer if the signal comes back during useage of the last memory
					rounds+= 1.
					mem = 0
			
			timeexit = time #this might be different that timesucc if we suceeded on the last memory
			memSucc = mem
			time2 = 0 # time after the first click
			

			if click == 0: #if the succ=1 is left from the previous failed purification attempt
				delta_t = time3exit - time3succ #time between the generation of the third click and stopping the previous purification attempt (see bottom)
				if succ == 1: #this is succ=1 left from the previous purification attempt
					if click2 ==0: #click2=0 corresponds to the case when in the previous attempts 3 copies were generated but the purification on the first 2 failed so the 3rd one was left (and we are still waiting for its herald)
						fid = fid3
					if click2 ==1: #this case corresponds to having generated two copies in the previous attempt, but while waiting for the herald of the second, the fidelity of the 1st one dropped below x and so the first one was discarded (we are still waiting for the herald of the second)
						fid = fid2
				if succ == 2: # this is the case if the 3 copies were generated in the previous attempt but while waiting for the heralds of 2 and 3, the fidelities of the first dropped below x so it was discarded, but we are still waiting for the heralds of 2 and 3 (so we have succ=2)
					fid = fid2
					fid2 = fid3
					


			while succ == 1 and time2 < d/c - Sg - (timeexit - timesucc) and click==1 and fid>x: #this is for the case when succ=1 was generated in this attempt, until the herald of the first successful attempts arrives
				time2 = time2 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt) #every time we try, the fidelity drops
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					succ += 1.
					time2succ = time2
				if mem == min(NoMem, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded
					time2 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c or more if the signal arrives while using the last memory
					rounds+= 1.
					mem = 0
				

			while succ == 1 and click==0 and time2 < d/c - Sg - delta_t and click2== 0 and fid>x: #this is for the case when succ=1 was left-over from the previous attempt (see above for what click2=0 means)
				time2 = time2 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt)
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					time2succ = time2
					succ += 1.
					click=1 #here we set it to one becuase the 2nd state would be generated the same way as when the first state is generated in this attempt
				if mem == min(NoMem, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded
					time2 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
					rounds+= 1.
					mem = 0
				

			while succ == 1 and click==0 and time2 < d/c - Sg - time3exit and click2 ==1 and fid>x: #this is for the case when succ=1 was left-over from the previous attempt (see above for what click2=1 means)
				time2 = time2 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt)
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					time2succ = time2
					succ += 1.
					click=1
				if mem == min(NoMem, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded
					time2 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
					rounds+= 1.
					mem = 0
				
			
			
			while succ == 1 and time2 >= d/c - Sg - (timeexit - timesucc) and click==1 and fid >= x: #this is for generating the second state after we already now that about the first success (first success in this attempt)
				time2 = time2 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt)
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					succ += 1.
					time2succ = time2
				if mem == min(NoMem - 1, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded, note that we have one less memory because the first is put aside to store
					time2 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
					rounds+= 1.
					mem = 0

			while succ == 1 and click==0 and time2 >= d/c - Sg - delta_t and click2 ==0 and fid >= x: #this is for generating the second state after we already now that about the first success (first success in the previous attempt, see above for def of click2=0)
				time2 = time2 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt)
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					succ += 1.
					click=1.
					time2succ = time2
				if mem == min(NoMem - 1, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded, note that we have one less memory because the first is put aside to store
					time2 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
					rounds+= 1.
					mem = 0
					

			while succ == 1 and click==0 and time2 >= d/c - Sg - time3exit and click2 ==1 and fid >= x:#this is for generating the second state after we already now that about the first success (first success in the previous attempt, see above for def of click2=1)
				time2 = time2 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt)
				if random.uniform(0, 1)< P * (1 - 0.5 * P):
					succ += 1.
					click=1.
					time2succ = time2
				if mem == min(NoMem - 1, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded, note that we have one less memory because the first is put aside to store
					time2 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
					rounds+= 1.
					mem = 0
					

			time2exit = time2 
			#at this moment we definitely have two clicks and now wait for the herald of the second
			time12 = time + time2 #total time we have spent already in this attempt
			memSucc = mem
			time3 = 0
			fid3 = 1 #we might succeed with the third copy while waiting for the herald of the second

			if click == 1: #this is if at least one of the successes was generated in this attempt
				fid2 = 1
				while fid>= x and time3 < d/c - Sg - (time2exit - time2succ): #we wait until receiving herald of the second success
					time3 = time3 + Eg + Sg
					mem += 1.
					fid =fid * numpy.exp(-1/mLt)
					fid2 = fid2 * numpy.exp(-1/mLt)
					if succ > 2:
						fid3 = fid3 * numpy.exp(-1/mLt)
					if random.uniform(0, 1)< P * (1 - 0.5 * P) and succ< 3: #we assume that getting 2 more successes while waiting for the herald of the second has negligible probability
						succ += 1.
						time3succ = time3
					if mem == min(NoMem - 1, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded
						time3 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time12 #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
						rounds+= 1.
						mem = 0
					

			#below is the case of waiting for the 2nd and 3rd heralds, if those two successed happened in the previous attempt. The third one should arrive at time3 = d/c -Sg - delta_t 
			while click ==0 and time3 < d/c - Sg - delta_t and fid>= x:
				time3 = time3 + Eg + Sg
				mem += 1.
				fid =fid * numpy.exp(-1/mLt)
				fid2 = fid2 * numpy.exp(-1/mLt)
				if succ > 2:
					fid3 = fid3 * numpy.exp(-1/mLt)
				if random.uniform(0, 1)< P * (1 - 0.5 * P) and succ< 3:
					succ += 1.
					time3succ = time3
				if mem == min(NoMem, numpy.ceil( d/(c*(Eg+Sg)) )): # if we have used all the memories but still haven't succeeded, note that we have one less memory because the first is put aside to store
					time3 = max( rounds * (d/c),rounds * mem * (Eg + Sg) ) - time12 #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c
					rounds+= 1.
					mem = 0
				

			time3exit = time3
			#now we have 2 heralded successes, possibly we have the third one, but don't know about it yet
			if fid >= x: #if the first copy still has good fidelity we do purification
				click2 = 0
				timefinal = timefinal + time12 + time3 + timepur + Sg  #our time is the time of the previous attempts + #time it takes to have 2 clicks in this attempt + time to receive the herald of the second + time to swap one of the successes to the electron spin + time of the controlled not gate and measurement
				succ -= 2.
				if random.uniform(0, 1) < 1./8.:
					succpur += 1.
					
			else:
				timefinal = timefinal + time12 + time3 #if the fidelity of the first copy is not good enough we need to discard it
				succ -= 1.
				click2 = 1.
				

			click=0 #mark that if we failed in this attempt and need to loop the whole attempt, then if we have some states left from this attempt we will use them in the next attempt

		timefinal = timefinal + d/c #if we suceeded in purification, we need to wait for the herald
		NoDec = min(numpy.ceil(d/c - Sg)/(Eg +Sg), NoMem-1) #this is the decoherence we get while waiting for the herald of teh successful purification
		decterm = fid * fid2 * numpy.exp( -NoDec/mLt )

		finalfid = (1+decterm)/2

		data[i,:] = [timefinal, finalfid]

	dataFixCut = numpy.mean(data, axis=0) #output: [duration, fidelity]
	
	return dataFixCut

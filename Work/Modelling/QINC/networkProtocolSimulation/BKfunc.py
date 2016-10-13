import random
import numpy

# Simulation of the multiplexed Barrett and Kok scheme, assuming that
# only one successful copy will be generated.
# Here the number of carbon memories is the input.
# That is we try until success and stop after receiving info about success.
# Hence we assume that we are not successful for the second time between the first successful run and the time it takes for the herald info to come back.
# This assumption can always be checked by comparing d/c (time of waiting for the herald during which we keep going) with "time" (duration of the whole protocol).
def BK(NoMem, n_iterations = 3000, pfc = 0.3,pout=0.3,d=50.0, Eg=1.,Sg=20.,mLt=400.):
	
	c=0.2

	if (NoMem-1) * (Eg + Sg) >= d/c: #If we have too many memories such that the heralding signal comes back before we would use all memories, then we will use only those first memories that can be all used before the signal comes back (or the signal can come back during the use of the last memory)
		NoMem = numpy.ceil(d/(c*(Eg+Sg)))
		print " Too many memories. We will only use ", NoMem,  " memories "

	
	Ploss = 10.**(-0.2*(d/20.)) #losses over the distance of d/2

	P = Ploss * pout * pfc
	data = numpy.zeros((n_iterations,2));
	for i in range(0, n_iterations):
		time = 0
		succ = 0
		mem = 0
		rounds = 1
		while succ == 0:
			time = time + Eg + Sg
			mem += 1.
			if random.uniform(0, 1)<0.5 * P**2:
				succ += 1.
			if mem == NoMem and succ==0: # if we have used all the memories but still haven't succeeded
				time = max( rounds * (d/c),rounds * NoMem * (Eg + Sg) ) #we wait till the signal for the first memory comes back, so that we finish a round lasting d/c, or if the signal arrived during using the last memory then the rounds finishes after we finish with that memory.
				rounds+= 1.
				mem = 0
		
		time = time + d/c - Sg #after the successful run we wait d/c to find out that it was successful, we substract Sg because d/c counts from the moment we finished Eg, before we do Sg
		# no of times we attempt to generate entanglement between the successful run and the heralding coming back. If NoMem is larger than (d/c)/(Eg + Sg) [that is NoMem = ceil((d/c)/(Eg + Sg))] then the first term is picked (which can be by 1 smaller than the no of memories, if NoMem is smaller than (d/c)/(Eg + Sg) then the second term is picked
		#basically first term corresponds to keep going for d/c - Sg time, and the second term, corresponds to trying as many times as we have memories and then waiting.
		NoDec = min(numpy.ceil(d/c - Sg)/(Eg +Sg), NoMem) #this is the number of trials that happen in d/c time after success. Depending on whehter we have more or less memories that d/c divided by the time of using one memory, we either try continuously in this window or if we have less memories, we try NoMem number of times and then wait.
		F = (1. + numpy.exp( -NoDec/mLt ))/2 #fidelity decays with this number of attempts after the successful run and before the herald info coming (when we stop)

		data[i,:] = [time, F]


	output = numpy.mean(data, axis=0)

	return [output[0], output[1]] #output: [duration, fidelity]

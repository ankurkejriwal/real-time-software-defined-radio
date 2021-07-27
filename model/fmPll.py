#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

# def fmPll(pllIn, freq, Fs, state, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01):

# 	# scale factors for proportional/integrator terms
# 	# these scale factors were derived assuming the following:
# 	# damping factor of 0.707 (1 over square root of 2)
# 	# there is no oscillator gain and no phase detector gain
# 	Cp = 2.666
# 	Ci = 3.555

# 	# gain for the proportional term
# 	Kp = (normBandwidth)*Cp
# 	# gain for the integrator term
# 	Ki = (normBandwidth*normBandwidth)*Ci

# 	# output array for the NCO
# 	ncoOut = np.empty(len(pllIn)+1)

# 	# initialize internal state
# 	# integrator = 0.0
# 	# phaseEst = 0.0
# 	# feedbackI = 1.0
# 	# feedbackQ = 0.0
# 	# ncoOut[0] = 1.0
# 	# trigOffset = 0
# 	# note: state saving will be needed for block processing
# 	ncoOut[0] = state['ncoLast']

# 	for k in range(len(pllIn)):

# 		# phase detector
# 		errorI = pllIn[k] * (+state['feedbackI'])  # complex conjugate of the
# 		errorQ = pllIn[k] * (-state['feedbackQ'])  # feedback complex exponential

# 		# four-quadrant arctangent discriminator for phase error detection
# 		errorD = math.atan2(errorQ, errorI)

# 		# loop filter
# 		state['integrator'] += Ki*errorD

# 		# update phase estimate
# 		state['phaseEst'] += Kp*errorD + state['integrator']
# 		# internal oscillator
# 		trigArg = 2*math.pi*(freq/Fs)*(state['trigOffset']+k+1) + state['phaseEst']
# 		state['feedbackI'] = math.cos(trigArg)
# 		state['feedbackQ'] = math.sin(trigArg)
# 		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)

# 	state['ncoLast'] = ncoOut[len(ncoOut)-1]
# 	state['trigOffset'] += len(pllIn) 
# 	# for stereo only the in-phase NCO component should be returned
# 	# for block processing you should also return the state
# 	return ncoOut, state
	# for RDS add also the quadrature NCO component to the output
import numpy as np
import math

def fmPll(pllIn, freq, Fs, recovery_state, ncoScale = 1.0, phaseAdjust = 0.0, \
    normBandwidth = 0.01):

    Cp = 2.666
    Ci = 3.555
    Kp = (normBandwidth)*Cp
    Ki = (normBandwidth*normBandwidth)*Ci

    ncoOut = np.empty(len(pllIn)+1)
    ncoOut_sin = np.empty(len(pllIn)+1)

    integrator = recovery_state[0]
    phaseEst = recovery_state[1]
    feedbackI = recovery_state[2]
    feedbackQ = recovery_state[3]
    ncoOut[0] = recovery_state[4]
    ncoOut_sin[0] = recovery_state[4]
    trigOffset = recovery_state[5]

    print(len(pllIn))

    for k in range(len(pllIn)):
        # phase detector
        errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
        errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential
        # four-quadrant arctangent discriminator for phase error detection
        errorD = math.atan2(errorQ, errorI)
        # loop filter
        integrator = integrator + Ki*errorD
        # update phase estimate
        phaseEst = phaseEst + Kp*errorD + integrator
        # internal oscillator
        trigArg = 2*math.pi*(freq/Fs)*(trigOffset+k+1) + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)
        ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
        ncoOut_sin[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)

    recovery_state[0] = integrator
    recovery_state[1] = phaseEst
    recovery_state[2] = feedbackI
    recovery_state[3] = feedbackQ
    recovery_state[4] = ncoOut[-1]
    recovery_state[5] = trigOffset + len(pllIn)

    return ncoOut, ncoOut_sin, recovery_state
	
if __name__ == "__main__":

	pass

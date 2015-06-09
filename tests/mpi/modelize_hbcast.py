#!/usr/bin/python

import sys
import matplotlib.pyplot as plt

from math import log

import hbcast_settings as st


class BcastModels:
    'Model based prediction of hbcast'

    def __init__(self, alpha = st.ALPHA, beta = st.BETA, msgSize = st.MSG_SIZE, msgSegment = st.SEGMENT_SIZE, alg = st.ALG):
        self.alg = alg
        self.msgSize = msgSize
        self.msgSegment = msgSegment
        self.numSegments = msgSize / float(msgSegment)
        self.alpha = alpha
        self.beta = beta

    def binomialBcast(self, communicatorSize):
        return log(communicatorSize, 2) * (self.alpha + self.msgSize * self.beta)

    def scRgAllgatherBcast(self, communicatorSize, numGroups):
	if (numGroups == 1 or numGroups == communicatorSize):
	    result = (log(communicatorSize, 2) + communicatorSize - 1) * self.alpha + 2 * self.msgSize * (1 - 1./communicatorSize) * self.beta
	else:
	    result = (log(communicatorSize, 2) + numGroups + communicatorSize/float(numGroups) - 1) * self.alpha + 2 * self.msgSize * (2 - 1./numGroups - float(numGroups)/communicatorSize) * self.beta
	    
	return result

    def scRdAllgatherBcast(self, communicatorSize, numGroups):
	if (numGroups == 1 or numGroups == communicatorSize):
            result = 2 * log(communicatorSize, 2) * self.alpha + 2 * self.msgSize * (1 - 1./communicatorSize) * self.beta
	else:
	    result = 2 * log(communicatorSize, 2) * self.alpha + 2 * self.msgSize * (2 - 1./numGroups - float(numGroups)/communicatorSize) * self.beta

	return result

    def linearBcast(self, communicatorSize):
        return (communicatorSize - 1) * (self.alpha + self.msgSize * self.beta)

    def pipelinedLinearBcast(self, communicatorSize):
        return (self.numSegments + communicatorSize - 2) * (self.alpha + self.msgSegment * self.beta)

    def binomialHbcast(self):
        #TODO
        pass

    def scRgAllgatherHbcast(self):
        #TODO
        pass

    def scRdAllgatherHBcast(self):
        #TODO
        pass

    def linearHbcast(self):
        #TODO
        pass

    def pipelinedLinearHbcast(self, communicatorSize, numGroups):
        if (numGroups == 1 or numGroups == communicatorSize):
            result = (self.numSegments + communicatorSize - 2) * (self.alpha + self.msgSegment * self.beta)
        else:
            result = (2*self.numSegments + (communicatorSize - numGroups*numGroups)/float(numGroups) - 4) * (self.alpha + self.msgSegment * self.beta)

        return result



    def pipelineBcast(self, communicatorSize, segmentSize):
        pass


    def mpichBcast(self, communicatorSize):
        if ( self.msgSize < st.SMALL_MSG_SIZE or communicatorSize <= st.MIN_COMM_SIZE ):
#	    print "run binomial %s " % (communicatorSize)
            return self.binomialBcast(communicatorSize)
	elif ( self.msgSize < st.LARGE_MSG_SIZE and self.isPowOf2(int(communicatorSize)) ):
#	    print "run scatter-rd-allgather %s " % (communicatorSize)
	    return bcastModels.scRdAllgatherBcast(communicatorSize, 1)
	else:
#	    print "run scatter-rg-allgather %s " % (communicatorSize)
	    return bcastModels.scRgAllgatherBcast(communicatorSize, 1)



    def mpichHbcast(self, communicatorSize, g):
#	print "mpichHbcast: p=%s, g=%s" % (p, g)
        return self.mpichBcast(communicatorSize/float(g)) + self.mpichBcast(g)



    def openMPIBcast(self, count, communicatorSize):
        smallMessageSize = 2048
        intermediateMessageSize = 370728
        a_p16  = 3.2118e-6
        b_p16  = 8.7936
        a_p64  = 2.3679e-6
        b_p64  = 1.1787
        a_p128 = 1.6134e-6
        b_p128 = 2.1102

        if ( self.msgSize < smallMessageSize or count <= 1 ):
            return self.binomialBcast(communicatorSize)
	elif ( self.msgSize < intermediateMessageSize ):
	    segSize = 1024
	    return self.splitBinaryBcast(communicatorSize, segSize)
	elif ( communicatorSize < 13 ):
	    segSize = 1024 << 3
	    return self.splitBinaryBcast(communicatorSize, segSize)
	elif ( communicatorSize < (a_p128 * self.msgSize + b_p128) ):
            segSize = 1024 << 7
            return self.pipelineBcast(communicatorSize, segSize)
	elif ( communicatorSize < (a_p64 * self.msgSize + b_p64) ):
	    segSize = 1024 << 6
	    return self.pipelineBcast(communicatorSize, segSize)
	elif ( communicatorSize < (a_p16 * self.msgSize + b_p16) ):
	    segSize = 1024 << 4
	    return self.pipelineBcast(communicatorSize, segSize)
	else:
	    segSize = 1024 << 3
	    return self.pipelineBcast(communicatorSize, segSize)




    def isPowOf2(self, number):
       return (number & (number-1) == 0) and (number != 0)


    def __str__(self):
        return "alg={:d}, msgSize={:f}, alpha={:f}, beta={:.12f}, segmentSize={:f}, numSegments={:f}".format(self.alg, self.msgSize, self.alpha, self.beta, self.msgSegment, self.numSegments)



#getcontext().prec = 6

if len(sys.argv) < 2:
    print "Usage: ./%s <p> <msg size>" % (sys.argv[0])
    sys.exit(-1)

numProc = int(sys.argv[1])
msg = float(sys.argv[2])

precision = st.PRECISION
minp=numProc
maxp=numProc
ming = 1
bcastModels = BcastModels(msgSize=msg)

binomialTimes = []
pipelinedTimes = []
scRgAllgatherTimes = []
scRdAllgatherTimes = []
mpichTimes = []
mpichHTimes = []
groups = []

p= maxp
for g in range(ming, p + 1):
    if (p % g == 0):
        groups.append(g)
	pipelinedTimes.append( round(bcastModels.pipelinedLinearBcast(p), precision) )
        binomialTimes.append( round(bcastModels.binomialBcast(p), precision) )
        scRgAllgatherTimes.append( round(bcastModels.scRgAllgatherBcast(p, g), precision) )
        scRdAllgatherTimes.append( round(bcastModels.scRdAllgatherBcast(p, g), precision) )
        mpichTimes.append( round(bcastModels.mpichBcast(p), precision) )
        mpichHTimes.append( round(bcastModels.mpichHbcast(p, g), precision) )
	

print "parameters    ", bcastModels
print "groups        ", groups
print "binomial      ", binomialTimes
print "pipelined     ", pipelinedTimes
print "scRgAllgather ", scRgAllgatherTimes
print "scRdAllgather ", scRdAllgatherTimes
print "mpichTimes    ", mpichTimes
print "mpichHTimes   ", mpichHTimes


#plotBinomial = plt.plot(groups, binomialTimes, color='green', label="binomial", marker="*")
#plotScRgAllgather = plt.plot(groups, scRgAllgatherTimes, color='black', label="Scatter-Rg-Allgather", marker="8")
plotScRdAllgather = plt.plot(groups, scRdAllgatherTimes, color='green', label="Scatter-Rd-Allgather", marker="^")
plotHMpich = plt.plot(groups, mpichHTimes, color='blue', label="hmpich", marker="+")
plotMpich =plt.plot(groups, mpichTimes, color='red', label='mpich', marker="s")
ax = plt.subplot(111)
ax.legend(loc='upper center',  bbox_to_anchor=(0., 1.02, 1., .102), fancybox=True, shadow=True, ncol=3)
plt.xscale('log', basex=2)
plt.xlabel('Number of groups')
plt.ylabel('Time in sec')
plt.title('m=%sKb, p=%s, alpha=%s, beta=%s' % (msg/1024., p, st.ALPHA, st.BETA), x=0.2, y=-0.1)
plt.show()





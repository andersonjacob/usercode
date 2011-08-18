maxDim = 11
maxDepth = 5
ebMaxPhi = 21
ebMaxEta = 86

def eta2ieta(eta):
    ieta = int(eta/0.087)
    if (ieta > 0):
        ieta += 1
    else:
        ieta -= 1
    return ieta

def phi2iphi(phi):
    iphi = int(phi/0.087) + 1
    return iphi

def HcalIndex(ieta, iphi, depth = 0):
    if (ieta < maxDim) and (iphi < maxDim) and (depth < maxDepth) and \
       (ieta > 0) and (iphi > 0) and (depth >= 0):
        if depth < 1:
            index = ieta*maxDim + iphi
        else:
            index = ieta*maxDim*maxDepth + iphi*maxDepth + depth
        ## print '({0},{1}) depth {2} => index: {3}'.format(ieta,iphi,
        ##                                                  depth,index)
        return index
    else:
        return -1

def EcalIndex(ieta, iphi):
    if (ieta < ebMaxEta) and (iphi < ebMaxPhi) and \
       (ieta > 0) and (iphi > 0):
        return ieta*ebMaxPhi + iphi
    else:
        return -1
    
def HcalEnergyAround(hits, ieta, iphi, depth = 0, radius = 1, calib = None):
    energy = 0.
    ## print "len hits:",len(hits)
    for e in range(ieta-radius,ieta+radius+1):
        for p in range(iphi-radius,iphi+radius+1):
            index = HcalIndex(e,p,depth)
            if (index >= 0) and (index < len(hits)):
                hitE = hits[index]
                if calib:
                    hitE *= calib[index]
                energy += hitE
    return energy

def EcalEnergyAround(hits, ieta, iphi, radius = 2, calib = None):
    energy = 0.
    for e in range(ieta-radius,ieta+radius+1):
        for p in range(iphi-radius,iphi+radius+1):
            index = EcalIndex(e,p)
            ## print '({0},{1}) => index: {2}'.format(e,p,index)
            if (index >= 0) and (index < len(hits)):
                hitE = hits[index]
                if calib:
                    hitE *= calib[index]
                energy += hitE
    return energy

def isInstrumented(det, ieta, iphi, depth, isHPD=False):
    if det.upper() == 'HB':
        return isInstrumentedHB(ieta, iphi, depth, isHPD)
    if det.upper() == 'HO':
        return isInstrumentedHO(ieta,iphi)
    return False

def isInstrumentedHB(ieta, iphi, depth, isHPD):
    if (iphi > 5) or (iphi < 2):
        return False
    if (iphi == 5) or (iphi == 2):
        isHPD = True
    if isHPD:
        if (ieta > 0) and (ieta < 15) and (depth < 2):
            return True
        if (ieta > 14) and (ieta < 17) and (depth < 3):
            return True
    elif (ieta > 5) and (ieta < 10) and (depth < 5):
        return True

    return False

def isInstrumentedHO(ieta, iphi):
    if (ieta < 5) and (ieta > 0) and \
       (iphi < 7) and (iphi > 2):
        return True
    elif (ieta < 11) and (ieta > 4) and \
         (iphi < 5) and (ieta > 1):
        return True
    return False

import cPickle as pickle

def loadCalibration(fileName):
    pkf = open(fileName, 'rb')
    newConst = pickle.load(pkf)
    pkf.close()
    return newConst

def storeCalibration(theConst, fileName):
    pkf = open(fileName, 'wb')
    pickle.dump(theConst, pkf)
    pkf.close()

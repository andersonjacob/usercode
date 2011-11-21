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

def digiIndex(ts, depth):
    if (depth < maxDepth) and (depth > 0) and \
           (ts < 10) and (ts >= 0):
        return ts*maxDepth + depth
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

def correctSaturation(pe, NP, xtalk=0.):
    from math import log
    if (pe >= NP):
        return 12.*NP
    if (pe < NP*0.05):
        return pe
    val = log(1.0 - float(pe)/NP)
    val *= -NP
    if (xtalk > 0) and (xtalk < 1.):
        val *= (1-xtalk)
    return val

def findBeamCenter(profile):
    return profile.GetBinCenter(profile.GetMaximumBin())

def findBeamCaloCoords(data):
    from ROOT import gDirectory
    data.Draw('maxEtaHB>>maxEtaHB(16,0.5,16.5)', 'triggerID==4', 'goff')
    temph = gDirectory.Get('maxEtaHB')
    ieta = findBeamCenter(temph)

    data.Draw('maxEtaEB>>maxEtaEB(50,0.5,50.5)', 'triggerID==4', 'goff')
    temph = gDirectory.Get('maxEtaEB')
    xtalieta = findBeamCenter(temph)

    data.Draw('maxPhiHB>>maxPhiHB(10,0.5,10.5)', 'triggerID==4', 'goff')
    temph = gDirectory.Get('maxPhiHB')
    iphi = findBeamCenter(temph)

    data.Draw('maxPhiEB>>maxPhiEB(21,0.5,21.5)', 'triggerID==4', 'goff')
    temph = gDirectory.Get('maxPhiEB')
    xtaliphi = findBeamCenter(temph)

    return (int(ieta), int(iphi), int(xtalieta), int(xtaliphi))

def loadConfig(filename):
    conf = {
        'HBmip1': 1.0,
        'HBmip2': 1.0,
        'HBmip3': 1.0,
        'HBmip4': 1.0,
        'HOmip': 1.0,
        'HBfCpe1': 1.0,
        'HBfCpe2': 1.0,
        'HBfCpe3': 1.0,
        'HBfCpe4': 1.0,
        'HOfCpe': 1.0,
        'beamE': 150.,
        'HBNP1': 0,
        'HBNP2': 0,
        'HBNP3': 0,
        'HBNP4': 0,
        'HONP': 0
        }
    if len(filename) > 0:
        for line in open(filename).readlines():
            tokens = line.rstrip().split('=')
            if (len(tokens) == 2):
                theKey = tokens[0].rstrip()
                theVal = tokens[1].rstrip()
                if theKey in conf:
                    conf[theKey] = type(conf[theKey])(theVal)
                else:
                    print 'unknown key',theKey,'with value',theVal
            else:
                print 'print invalid line:',line.rstrip()

    return conf

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

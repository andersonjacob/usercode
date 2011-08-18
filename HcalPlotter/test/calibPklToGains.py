from tbRoutines import *

from optparse import OptionParser

parser = OptionParser()

(opts,args) = parser.parse_args()

calibConst = loadCalibration(args[0])

calibFile = open(args[0].replace('pkl', 'txt'), 'w')

calibFile.write('#U ADC\n')
calibFile.write('# eta  phi  dep  det  cap1    cap2    cap3    cap4     HcalDetId\n')
theLine = '  {0}    {1}   {2}    {3}   {4:0.3f}   {4:0.3f}   {4:0.3f}   {4:0.3f}    {5:x}\n'

for eta in range(1, 17):
    for phi in range(1, 8):
        for depth in range(1, maxDepth):
            if isInstrumented('HB', eta, phi, depth):
                if (eta < 15):
                    encodePhi = depth*10+phi
                    encodeDepth = 1
                else:
                    encodePhi = 10+phi
                    encodeDepth = depth
                index = HcalIndex(eta, phi, depth)
                if (index > 0):
                    calibFile.write(theLine.format(eta,encodePhi,encodeDepth,
                                                   'HB',calibConst[index],0))
                else:
                    calibFile.write(theLine.format(eta,encodePhi,encodeDepth,
                                                   'HB',1.0,0))

calibFile.write('# HO\n')
for eta in range(1, 17):
    for phi in range(1, 8):
        if isInstrumented('HO', eta, phi, 4):
            calibFile.write(theLine.format(eta,phi,4,'HO',1.0,0))

calibFile.close()
                

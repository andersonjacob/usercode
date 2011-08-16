#! /usr/bin/python

import os

def createFileName (eta, bxints, NP, type, test, test2, rTime1, rPct1, rTime2):
    return "HistEta%dInts%dPx%dType%dTest%dTest2%dT1t%dTL%3.2fRt%d" % (eta, bxints, NP, type, test, test2, rTime1, rPct1, rTime2)

def RunHistory(eta, bxints, NP, type, test, rTime1, rPct1, rTime2, pde, past,
               test2=0):
    logfile = createFileName(eta, bxints, NP, type, test, test2, rTime1,
                             rPct1, rTime2)
    cmd = 'root -l -b -q \'SiPMHistoryStudy.cc+(100,' + str(eta) + \
          ',' + str(bxints) + ',' + str(NP) + ',' + str(type) + ',' + \
          str(past) + ',' + str(rTime1) + ',' + str(rPct1) + ',' + \
          str(rTime2) + ',' + logfile + '.asc,' + str(pde) + \
          ',0,\"peSpectrumByLayer7EtaBins.root\",\"LayerPes' + str(test) + \
          'GeV.txt\",' + str(test2) +')\' 1> ' + logfile + '.log'

    print cmd
    os.system(cmd)

def PlotHistory(eta, bxints, NP, type, test, rTime1, rPct1, rTime2,
                test2=0):
    fname = createFileName(eta, bxints, NP, type, test, test2, rTime1,
                           rPct1, rTime2) + '.asc'
    tpe = "Type%dPx%dEta%dT1t%dTL%3.2fRt%d" % (type, NP, eta, rTime1, rPct1, rTime2)
    cmd = 'root -l -b -q \'plotSiPMHistory.cc+(\"' + fname + \
          '\",\"plots/' + tpe + '\",' + str(bxints) + ',' + str(test) + ',' + \
          str(test2) + ',\"tables/' + tpe + 'PxTab.tex\",\"tables/' + tpe + \
          'RespTab.tex\",\"tables/' + tpe + 'OccTab.tex\")\''
    print cmd
    os.system("mkdir -p plots/" + tpe)
    os.system(cmd)
    
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--etaBin", type="int", dest="etaBin", metavar="ETABIN",
                      default=1, help="which eta bin between 1 and 7");
    parser.add_option("--bxpu", type="int", dest="intsperx", metavar="PILEUP",
                      default=1, help="average number of interactions per bunch crossing");
    parser.add_option("--edu", action="store_true", dest="edu", default=False);
    parser.add_option("--odu", action="store_true", dest="odu", default=False);
    parser.add_option("--NP", type="int", dest="NP", metavar="PIXELDENSITY",
                      default=15000, help="pixels per mm^2");
    parser.add_option("--PDE", type="float", dest="pde", metavar="RELPDE",
                      default=1.0, help="photon detection efficiency relative to 400/mm^2 Hamamatsu");
    parser.add_option("--past", type="int", dest="past", metavar="PASTLEN",
                      default=1000000, help="how man ns back in time to go");
    parser.add_option("--tTime", type="float", dest="rTime1",
                      metavar="RECOVERYTRANSITIONTIME", default=160.0,
                      help="time in ns of the elbow in recovery")
    parser.add_option("--tLevel", type="float", dest="rPct1",
                      metavar="RECOVERYTRANSITIONLEVEL", default=0.47,
                      help="recharge level at the elbow between 0 and 1")
    parser.add_option("--tFull", type="float", dest="rTime2",
                      metavar="FULLRECOVERYTIME", default=1E6,
                      help="time in ns of full recovery");
    parser.add_option("--test", type="int", dest="test", metavar="TESTGEV",
                      default=500, help="which test file to use (1000 or 500)")
    parser.add_option("--test2", type="int", dest="test2", metavar="TESTGEV2",
                      default=1000, help="which second test file to use (1000 or 500), 0 for none")
    

    (options, args) = parser.parse_args()

    if (options.edu is options.odu):
        parser.error("You must have only either an edu or an odu.");

    print options

    SiPMType = 0
    if options.edu:
        SiPMType = 1
    RunHistory(options.etaBin, options.intsperx, options.NP, SiPMType,
               options.test, options.rTime1, options.rPct1, options.rTime2,
               options.pde, options.past, options.test2)
    PlotHistory(options.etaBin, options.intsperx, options.NP, SiPMType,
                options.test, options.rTime1, options.rPct1, options.rTime2,
                options.test2)
        

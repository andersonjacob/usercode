from ROOT import TFitter, gROOT
from array import array

gROOT.ProcessLine('.L truncRMS.cc+')

from ROOT import setData, hookupMinuit

def runMinimization(calibData, MPV, errLevel = 0.001):
    for param in calibData.keys():
        rowData = array('d', calibData[param])
        print param,'len:',len(calibData[param]) #, param[0:2], param[:2]
        setData(len(calibData[param]), rowData, param)

    minner = TFitter(len(calibData.keys())*2+1)

    pl = array('d', [1.])
    errdef = array('d', [errLevel])
    minner.ExecuteCommand('SET PRINT', pl, 1)
    minner.ExecuteCommand('SET ERR', errdef, 1)

    hookupMinuit(minner)

    parN = 0
    minner.SetParameter(parN, "MPV", MPV, 10., 0., 500.)
    minner.FixParameter(0)
    parN += 1
    for param in calibData.keys():
        initVal = 1.0
        minner.SetParameter(parN, param, initVal, initVal/10., 0., 10.)
        parN += 1

    arglist = array('d', [len(calibData.keys())*500., 0.1])
    ## minner.ExecuteCommand('SIMPLEX', arglist, 2)
    arglist[1] = 1.
    minner.ExecuteCommand('MINIMIZE', arglist, 2)
    ## minner.ExecuteCommand('HESSE', arglist, 1)
    return minner

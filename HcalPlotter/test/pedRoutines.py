from ROOT import RooDataSet, RooRealVar, RooGaussian, RooArgSet, RooFit,\
     gPad, RooAddPdf, RooWorkspace, RooFormulaVar, RooArgList, RooMsgService,\
     RooBifurGauss

RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

def fillDataSet(data, x, N, dsName = 'ds'):
    cols = RooArgSet(x)
    ds = RooDataSet(dsName, dsName, cols)
    #ds.Print()
    print 'length data:', N
    for datum in range(0,N):
        if (data[datum] < x.getMax()) and (data[datum] > x.getMin()):
            x.setVal(data[datum])
            ds.add(cols)
    ds.Print()
    return ds

def fitPed(hist, ws, name='x'):
    maxBin = hist.GetMaximumBin()
    x = ws.var(name)
    #rds = ds.reduce('{1}<{0:0.2f}'.format(hist.GetBinLowEdge(maxBin+2),name))
    #rds.Print()
    x.setRange('ped_fit', x.getMin(), hist.GetBinLowEdge(maxBin+3))
    pedMean = RooRealVar('pedMean', 'mean_{ped}', hist.GetBinCenter(maxBin),
                         x.getMin(), x.getMax())
    pedMean.Print()
    pedWidth = RooRealVar('pedWidth', 'sigma_{ped}', 1., 0., 10.)
    pedWidth.Print()
    ped = RooGaussian('ped', 'ped', x, pedMean, pedWidth)

    pedMean.setConstant(False)
    ped.fitTo(ws.data('ds'), RooFit.Minos(False), RooFit.Range('ped_fit'),
              RooFit.PrintLevel(0))

    getattr(ws, 'import')(ped)
    ## xf = x.frame(x.getMin('ped_fit'), x.getMax('ped_fit'),
    ##              int(x.getMax('ped_fit')-x.getMin('ped_fit')))
    ## ws.data('ds').plotOn(xf)
    ## ped.plotOn(xf)
    ## xf.Draw()
    ## gPad.Update()
    ## gPad.WaitPrimitive()

def findOnePe(hist, ws, name='x', Npe = 1):
    fitPed(hist, ws, name)
    x = ws.var(name)

    ped = ws.pdf('ped')
    pedWidth = ws.var('pedWidth')

    pdfs = RooArgList(ped)
    pdfList = []

    fped = RooRealVar('fped', 'f_{ped}', 0.8, 0., 1.)
    fractions = RooArgList(fped)
    fList = []
    peList = []

    peMean = RooRealVar('peMean', 'mean_{pe}', 6., 0., 20.)
    peWidth = RooRealVar('peWidth', 'width_{pe}', pedWidth.getVal(), 0., 10.)

    for i in range(0, Npe):
        pem = RooFormulaVar('pem{0}'.format(i+1), '@0+{0}*@1'.format(i+1),
                            RooArgList(ws.var('pedMean'), peMean))
        peList.append(pem)
        npepdf = RooGaussian('pe{0}pdf'.format(i+1), 'pe{0}pdf'.format(i+1),
                             x, pem, pedWidth)
        pdfs.add(npepdf)
        pdfList.append(npepdf)

        fnpe = RooRealVar('fpe{0}'.format(i+1), 'fpe{0}'.format(i+1),
                          0.5, -0.1, 1.0)
        fractions.add(fnpe)
        fList.append(fnpe)

    #bgMean = RooRealVar("bgMean", "bgMean", 6.0, x.getMin(), x.getMax())
    bgScale = RooRealVar("bgScale", "bgScale", 0.5, -1.0, Npe + 1.0)
    bgMean = RooFormulaVar("bgMean", "@1+@0*@2",
                           RooArgList(peMean, ws.var('pedMean'), bgScale))
    bgWidthL = RooRealVar("bgWidthL", "bgWidthL", pedWidth.getVal()*2,
                          0., 25.)
    bgWidthR = RooRealVar("bgWidthR", "bgWidthR", pedWidth.getVal()*7,
                          0., 25.)

    bgGauss = RooBifurGauss("bgGauss", "bgGauss", x, bgMean,
                            bgWidthR, bgWidthR)

    if (Npe > 1):
        pdfs.add(bgGauss)
    else:
        fractions.remove(fractions.at(fractions.getSize()-1))

##     pem = RooFormulaVar('pem', '@0+@1', RooArgList(peMean, ws.var('pedMean')))
##     firstPe = RooGaussian('firstPe', 'firstPe', x, pem, peWidth)

##     pdfs.Print("v")
##     fractions.Print("v")
    pedPlusOne = RooAddPdf('pedPlusOne', 'pedPlusOne', pdfs, fractions, True)

    ## pedWidth = ped.GetParameter(2)
    ## pedMean = ped.GetParameter(1)
    ## pedA = ped.GetParameter(0)
    
    secondMax = hist.GetMaximumBin() + 1
    goingDown = True
    maxVal = hist.GetBinContent(secondMax)
    foundMax = False
    while (not foundMax) and (secondMax < hist.GetNbinsX()):
        tmpVal = hist.GetBinContent(secondMax+1)
        if (tmpVal < maxVal):
            if not goingDown:
                foundMax = True
            else:
                goingDown = True
                maxVal = tmpVal
                secondMax += 1
        elif (tmpVal > maxVal):
            goingDown = False
            maxVal = tmpVal
            secondMax += 1
        else:
            maxVal = tmpVal
            secondMax += 1

    secondMaxx = hist.GetBinCenter(secondMax)
    print 'found 2nd maximum in bin',secondMax,'value',secondMaxx

##     peMean.setVal(secondMaxx)
##     bgMean.setVal(secondMaxx*0.6)
    x.setRange('pedPlus_fit', x.getMin(), ws.var('pedMean').getVal()+pedWidth.getVal()*6.*(Npe+0))

    pedPlusOne.fitTo(ws.data('ds'), RooFit.Minos(False),
                     RooFit.Range('pedPlus_fit'),
                     RooFit.PrintLevel(1))

    getattr(ws, 'import')(pedPlusOne)

if __name__ == '__main__':
    import sys, os
    sys.path.append(os.environ['HOME']+'/pyroot')
    del sys
    del os
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--phi', type='int', default=2, dest='phi',
                      help='override table iphi position')
    parser.add_option('--eta', type='int', default=5, dest='eta',
                      help='override table ieta position')
    parser.add_option('--npe', type='int', default=1, dest='npe',
                      help='Number of pe to try to fit')
    (opts, args) = parser.parse_args()

    import root_logon

    from ROOT import TFile, TTree, gDirectory, TMath, kRed, kDashed

    inFile = TFile(args[0])
    dataTree = inFile.Get("plotanal/dataTree")

    HBDigis = int(dataTree.GetMaximum('NHBdigis'))
    HODigis = int(dataTree.GetMaximum('NHOdigis'))

    pedCut = '(triggerID==1)&&(NHOdigis=={0})'.format(HODigis)
    HOTower = 'HOE[{0}][{1}]'.format(opts.eta, opts.phi)
    dataTree.Draw(HOTower, pedCut, 'goff')

    print 'selected rows:',dataTree.GetSelectedRows(),dataTree.GetV1()
    
    minPed = TMath.MinElement(dataTree.GetSelectedRows(), dataTree.GetV1())
    maxPed = TMath.MaxElement(dataTree.GetSelectedRows(), dataTree.GetV1())
    #print 'minPed:',minPed,'maxPed:',maxPed
    minPed = int(minPed) - 2.5
    maxPed = int(maxPed) + 2.5
    print 'minPed:',minPed,'maxPed:',maxPed
    while int((maxPed-minPed)/2.)*2. < (maxPed-minPed):
        maxPed += 1.
    print 'minPed:',minPed,'maxPed:',maxPed
    dataTree.Draw('{0}>>pedhist({1},{2:0.1f},{3:0.1f}'.format(HOTower,
                                                              int(maxPed-minPed),
                                                              minPed,maxPed),
                  pedCut, 'goff')

    pedhist = gDirectory.Get('pedhist')
    ws = RooWorkspace('ws')
    x = RooRealVar('x', 'energy', minPed, maxPed, 'fC')
    x.Print()

    ds = fillDataSet(dataTree.GetV1(), x, dataTree.GetSelectedRows())
    getattr(ws, 'import')(ds)
    findOnePe(pedhist, ws, Npe=opts.npe)

    pedPlusOne = ws.pdf('pedPlusOne')
    peMean = ws.var('peMean')
    pedMean = ws.var('pedMean')
    ## pedhist.Draw()
    ## onePeF.Draw('same')

    xf = x.frame(x.getMin(), x.getMax(), int(x.getMax()-x.getMin()))
    ds.plotOn(xf)
    pedPlusOne.plotOn(xf) #,
                      #RooFit.Range('pedPlus_fit'))
    pedPlusOne.plotOn(xf, #RooFit.Range('pedPlus_fit'),
                      RooFit.Components("bg*"),
                      RooFit.LineColor(kRed),
                      RooFit.LineStyle(kDashed))
    xf.Draw()
    gPad.Update()

    print '\nsingle fC/pe for HO ({0}, {1}):'.format(opts.eta, opts.phi),\
          '{0:0.3f} +/- {1:0.3f}'.format(peMean.getVal(), peMean.getError())

from ROOT import RooDataSet, RooRealVar, RooGaussian, RooArgSet, RooFit,\
     gPad, RooAddPdf, RooWorkspace, RooFormulaVar, RooArgList

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

def findOnePe(hist, ws, name='x'):
    fitPed(hist, ws, name)
    x = ws.var(name)

    ped = ws.pdf('ped')
    pedWidth = ws.var('pedWidth')

    peMean = RooRealVar('peMean', 'mean_{pe}', 4., 0., 20.)
    pem = RooFormulaVar('pem', '@0+@1', RooArgList(peMean, ws.var('pedMean')))
    peWidth = RooRealVar('peWidth', 'width_{pe}', pedWidth.getVal(), 0., 10.)
    firstPe = RooGaussian('firstPe', 'firstPe', x, pem, pedWidth)

    fped = RooRealVar('fped', 'f_{ped}', 0.8, 0., 1.)

    pedPlusOne = RooAddPdf('pedPlusOne', 'pedPlusOne', ped, firstPe, fped)

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

    x.setRange('pedPlus_fit', x.getMin(), ws.var('pedMean').getVal()+pedWidth.getVal()*6.)

    pedPlusOne.fitTo(ws.data('ds'), RooFit.Minos(False),
                     RooFit.Range('pedPlus_fit'), RooFit.PrintLevel(0))

    getattr(ws, 'import')(pedPlusOne)
    ## xf = x.frame(x.getMin('pedPlus_fit'), x.getMax('pedPlus_fit'),
    ##              int(x.getMax('pedPlus_fit')-x.getMin('pedPlus_fit')))
    ## ds.plotOn(xf)
    ## pedPlusOne.plotOn(xf)
    ## xf.Draw()
    ## gPad.Update()
    ## gPad.WaitPrimitive()

    #return pedPlusOne, ped, pedMean, peMean

    ## pedPlusOne = TF1("pedPlusOne", "gaus(0)+gaus(3)", hist.GetBinLowEdge(1),
    ##                  secondMaxx+3.0*pedWidth)
    ## for p in range(0,3):
    ##     pedPlusOne.SetParameter(p, ped.GetParameter(p));

    ## pedPlusOne.SetParameter(3, maxVal)
    ## pedPlusOne.SetParameter(4, secondMaxx)
    ## pedPlusOne.SetParameter(5, pedWidth)

    ## #pedPlusOne.Print()

    ## hist.Fit(pedPlusOne, 'LR0E')
    ## return pedPlusOne

if __name__ == '__main__':
    import sys, os
    sys.path.append(os.environ['HOME']+'/pyroot')
    del sys
    del os
    import root_logon

    from ROOT import TFile, TTree, gDirectory, TMath
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--phi', type='int', default=2, dest='phi',
                      help='override table iphi position')
    parser.add_option('--eta', type='int', default=5, dest='eta',
                      help='override table ieta position')
    (opts, args) = parser.parse_args()

    inFile = TFile(args[0])
    dataTree = inFile.Get("plotanal/dataTree")

    pedCut = '(triggerID==1)&&(NHOdigis==33)'
    HOTower = 'HOE[{0}][{1}]'.format(opts.eta, opts.phi)
    dataTree.Draw(HOTower, pedCut, 'goff')

    print 'selected rows:',dataTree.GetSelectedRows(),dataTree.GetV1()
    
    minPed = TMath.MinElement(dataTree.GetSelectedRows(), dataTree.GetV1())
    maxPed = TMath.MaxElement(dataTree.GetSelectedRows(), dataTree.GetV1())
    #print 'minPed:',minPed,'maxPed:',maxPed
    minPed = int(minPed) - 2.5
    maxPed = int(maxPed) + 2.5
    if (int(maxPed-minPed)/2)*2 < (maxPed-minPed):
        maxPed += 1
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
    findOnePe(pedhist, ws)

    pedPlusOne = ws.pdf('pedPlusOne')
    peMean = ws.var('peMean')
    pedMean = ws.var('pedMean')
    ## pedhist.Draw()
    ## onePeF.Draw('same')

    xf = x.frame(x.getMin(), x.getMax(), int(x.getMax()-x.getMin()))
    ds.plotOn(xf)
    pedPlusOne.plotOn(xf, RooFit.Range('pedPlus_fit'))
    xf.Draw()
    gPad.Update()

    print '\nsingle fC/pe for HO ({0}, {1}):'.format(opts.eta, opts.phi),\
          '{0:0.2f}'.format(peMean.getVal()-pedMean.getVal())

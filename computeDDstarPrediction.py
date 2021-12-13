'''
Script for the computation of the expected D-D* correlation function
'''

import os
import sys
import argparse
import pandas as pd
import numpy as np
import yaml
from ROOT import TFile, TGaxis, TCanvas, TF1, TGraph, TH1F, TSpline3, TLegend, TLatex, TLine # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kGreen, kOrange, kRed, kBlack, kGray, kFullCircle, kOpenCircle, kFullDiamond, gStyle, kRainBow # pylint: disable=import-error,no-name-in-module

gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadLeftMargin(0.12)
gStyle.SetPadRightMargin(0.035)
gStyle.SetTitleOffset(1.2, 'yz')
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(kRainBow)
TGaxis.SetMaxDigits(3)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', metavar='text', default='config.yml')
parser.add_argument('--yMax', type=float, default=8.)
parser.add_argument('--xMax', type=float, default=0.5)
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlcfgFileName:
    cfg = yaml.load(ymlcfgFileName, yaml.FullLoader)

inDirName = cfg['inputdir']
doAccStudy = cfg['rapidity']['plot']
Dspecie = cfg['species']['D']
Dstarspecie = cfg['species']['Dstar']

if Dspecie == 411:
    Dtitle = 'D^{#plus}'
    brD = 0.0938
elif Dspecie == -411:
    Dtitle = 'D^{#minus}'
    brD = 0.0938
elif Dspecie == 421:
    Dtitle = 'D^{0}'
    brD = 0.03951
elif Dspecie == -421:
    Dtitle = '#bar{D}^{0}'
    brD = 0.03951
else:
    print(f'ERROR: D specie {Dspecie} not supported! Exit')
    sys.exit()

if Dstarspecie == 413:
    Dstartitle = 'D*^{#plus}'
    brDstar = 0.677 * 0.03951
elif Dstarspecie == -413:
    Dstartitle = 'D*^{#minus}'
    brDstar = 0.677 * 0.03951
elif Dstarspecie == 423:
    Dstartitle = 'D*^{0}'
    brDstar = 0.353 * 0.03951
elif Dstarspecie == -423:
    Dstartitle = '#bar{D}*^{0}'
    brDstar = 0.353 * 0.03951
else:
    print(f'ERROR: D* specie {Dstarspecie} not supported! Exit')
    sys.exit()

inFileNames = os.listdir(inDirName)
nEvents = cfg['eventsperfile'] * len(inFileNames)

print(f'\nAnalysing simulation outputs with {len(inFileNames)} files for a total of {nEvents} events\n')

hSEDistrVsPt, hMEDistrVsPt, hSEDistrVsY = None, None, None
for inFileName in inFileNames:
    inFile = TFile.Open(os.path.join(inDirName, inFileName))
    hTmpSE = inFile.Get(f'hPairSE_{Dstarspecie}_{Dspecie}')
    hTmpME = inFile.Get(f'hPairME_{Dstarspecie}_{Dspecie}')
    hTmpSE.SetDirectory(0)
    hTmpME.SetDirectory(0)
    if doAccStudy:
        hTmpVsY = inFile.Get(f'hPairVsY_{Dstarspecie}_{Dspecie}')
    if hSEDistrVsPt:
        hSEDistrVsPt.Add(hTmpSE)
        hMEDistrVsPt.Add(hTmpME)
        if doAccStudy:
            hSEDistrVsY.Add(hTmpVsY)
    else:
        hSEDistrVsPt = hTmpSE
        hMEDistrVsPt = hTmpME
        if doAccStudy:
            hSEDistrVsY = hTmpVsY
    inFile.Close()

BRfactor = brD * brDstar
lumiScaleFactorPP = 6.e14 / nEvents
lumiScaleFactorPbPb = 12.5e9 * (1799.9**2 + 1405.9**2) / nEvents

hSEDistr = hSEDistrVsPt.ProjectionZ()
hMEDistr = hMEDistrVsPt.ProjectionZ()
hSEDistr.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')
hSEDistr.GetXaxis().SetTitleSize(0.045)
hSEDistr.GetXaxis().SetLabelSize(0.045)
hMEDistr.GetXaxis().SetTitleSize(0.045)
hMEDistr.GetXaxis().SetLabelSize(0.045)
hSEDistr.GetYaxis().SetTitle('#it{N}_{pairs}^{same}')
hMEDistr.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')
hMEDistr.GetYaxis().SetTitle('1/#it{N} #it{N}_{pairs}^{mixed}')

hPtDvsKstar = hSEDistrVsPt.Project3D('xz')
hPtDstarvsKstar = hSEDistrVsPt.Project3D('yz')
hPtDvsKstar.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')
hPtDvsKstar.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')
hPtDstarvsKstar.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')
hPtDstarvsKstar.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')

inFileDmes = TFile.Open(cfg['Dmesonperf']['inputfile'])
ptMinsD = cfg['Dmesonperf']['ptmins']
ptMaxsD = cfg['Dmesonperf']['ptmaxs']
ptLimits = ptMinsD.copy()
ptLimits.append(ptMaxsD[-1])
hSignalOverBkg = TH1F('hSignalOverBkg', '', len(ptMaxsD), np.array(ptLimits, 'd'))
fGaus = TF1('fGaus', 'gaus', 1.8, 2.)
for iPtD, (ptMinD, ptMaxD) in enumerate(zip(ptMinsD, ptMaxsD)):
    hSignal = inFileDmes.Get(f'invmass_signal_rap0_pt{iPtD}')
    hBkg = inFileDmes.Get(f'invmass_background_rap0_pt{iPtD}')
    hSignal.Fit('fGaus', 'Q0')
    mean = fGaus.GetParameter(1)
    sigma = fGaus.GetParameter(2)
    minMassBin = hSignal.GetXaxis().FindBin(mean-3*sigma)
    maxMassBin = hSignal.GetXaxis().FindBin(mean+3*sigma)
    signal = hSignal.Integral(minMassBin, maxMassBin)
    bkg = hBkg.Integral(minMassBin, maxMassBin)
    hSignalOverBkg.SetBinContent(iPtD+1, signal/bkg)

    if iPtD == 0:
        hTot = inFileDmes.Get(f'invmass_signalbackground_rap0_pt{iPtD}')
        hTot.SetLineColor(kBlack)
        hTot.SetLineWidth(2)
        cMass = TCanvas('cMass', '', 800, 800)
        cMass.DrawFrame(1.72, 0., 2.05, 4000, ';#it{M}_{K#pi} (GeV/#it{c}) ;entries')
        hTot.Draw('same')
        cMass.Modified()
        cMass.Update()

hSignalOverBkgDVsKstar = hSEDistr.Clone('hSignalOverBkgDVsKstar')
hSignalOverBkgDstarVsKstar = hSEDistr.Clone('hSignalOverBkgDstarVsKstar')
nKstarBins = hPtDvsKstar.GetXaxis().GetNbins()
nPtBins = hPtDvsKstar.GetYaxis().GetNbins()
for iKstarBin in range(1, nKstarBins+1):
    SoverBKstarBin = 0
    for iPtBin in range(1, nPtBins+1):
        ptCent = hPtDvsKstar.GetYaxis().GetBinCenter(iPtBin)
        ptBinSoB = hSignalOverBkg.GetXaxis().FindBin(ptCent)
        SoverBKstarBin += hSignalOverBkg.GetBinContent(ptBinSoB)
    SoverBKstarBin /= nPtBins
    hSignalOverBkgDVsKstar.SetBinContent(iKstarBin, SoverBKstarBin)
for iKstarBin in range(1, nKstarBins+1):
    SoverBKstarBin = 0
    for iPtBin in range(1, nPtBins+1):
        ptCent = hPtDstarvsKstar.GetYaxis().GetBinCenter(iPtBin)
        ptBinSoB = hSignalOverBkg.GetXaxis().FindBin(ptCent)
        SoverBKstarBin += hSignalOverBkg.GetBinContent(ptBinSoB)
    SoverBKstarBin /= nPtBins
    hSignalOverBkgDstarVsKstar.SetBinContent(iKstarBin, SoverBKstarBin)

cPtVsKstar = TCanvas('cPtVsKstar', '', 1200, 600)
cPtVsKstar.Divide(2, 1)
cPtVsKstar.cd(1)
hPtDvsKstar.GetXaxis().SetTitleSize(0.045)
hPtDvsKstar.GetXaxis().SetLabelSize(0.045)
hPtDvsKstar.GetYaxis().SetTitleSize(0.045)
hPtDvsKstar.GetYaxis().SetLabelSize(0.045)
hPtDvsKstar.Draw('colz')
cPtVsKstar.cd(2)
hPtDstarvsKstar.Draw('colz')
hPtDstarvsKstar.GetXaxis().SetTitleSize(0.045)
hPtDstarvsKstar.GetXaxis().SetLabelSize(0.045)
hPtDstarvsKstar.GetYaxis().SetTitleSize(0.045)
hPtDstarvsKstar.GetYaxis().SetLabelSize(0.045)
cPtVsKstar.Modified()
cPtVsKstar.Update()

hSEDistr.Scale(lumiScaleFactorPP * BRfactor)
hMEDistr.Scale(lumiScaleFactorPP * BRfactor)
hSEDistr.Rebin(5)
hMEDistr.Rebin(5)
hSEDistr.SetLineColor(kBlack)
hMEDistr.SetLineColor(kBlack)
hSEDistr.SetLineWidth(2)
hMEDistr.SetLineWidth(2)
hSEDistr.SetMarkerColor(kBlack)
hMEDistr.SetMarkerColor(kBlack)
hSEDistr.SetMarkerStyle(kFullCircle)
hMEDistr.SetMarkerStyle(kFullCircle)
hSignalOverBkgDVsKstar.Rebin(5)
hSignalOverBkgDstarVsKstar.Rebin(5)
hSignalOverBkgDVsKstar.Scale(1./5)
hSignalOverBkgDstarVsKstar.Scale(1./5)

hCFPythia = hSEDistr.Clone('hCFPythia')
hCFPythia.GetYaxis().SetDecimals()
hCFPythia.GetYaxis().SetTitle(f'#it{{C}}_{{{Dtitle}{Dstartitle}}}')
hCFPythia.Divide(hMEDistr)
fPol1 = TF1('fPol1', 'pol1', 0., 2.)
fPol1.SetLineWidth(2)
fPol1.SetLineColor(kAzure+4)
hCFPythia.Fit('fPol1', '0', '', 0.2, 2.)
isFlat = False
if abs(fPol1.GetParameter(1)) < 1.e-4: # small deviation from 1, negligible
    fPol0 = TF1('fPol0', 'pol0', 0., 2.)
    fPol0.SetLineWidth(2)
    fPol0.SetLineColor(kAzure+4)
    hCFPythia.Fit('fPol0', '0', '', 0.2, 2.)
    scaleFact = fPol0.GetParameter(0)
    isFlat = True # no jet background
    fPol0.SetParameter(0, fPol0.GetParameter(0)/scaleFact)
    hCFPythia.Scale(1./scaleFact)
else:
    scaleFact = fPol1.Eval(2)
    hCFPythia.Scale(1./scaleFact)
    fPol2 = TF1('fPol2', 'pol2', 0., 2.)
    fPol2.SetLineWidth(2)
    fPol2.SetLineColor(kAzure+4)
    hCFPythia.Fit('fPol2', '0', '', 0.2, 2.)
hMEDistr.Scale(scaleFact)

lineAtOne = TLine(0., 1., 2., 1.)
lineAtOne.SetLineWidth(2)
lineAtOne.SetLineStyle(9)
lineAtOne.SetLineColor(kGray+1)

cDistrPythia = TCanvas('cDistrPythia', '', 1250, 500)
cDistrPythia.Divide(3, 1)
cDistrPythia.cd(1).SetLogy()
hSEDistr.Draw('e')
cDistrPythia.cd(2).SetLogy()
hMEDistr.Draw('e')
cDistrPythia.cd(3)
if hCFPythia.GetMaximum() > 8.:
    hCFPythia.GetYaxis().SetRangeUser(0., 8.)
hCFPythia.Draw('e')
lineAtOne.Draw()
if not isFlat:
    fPol2.Draw('same')
cDistrPythia.Modified()
cDistrPythia.Update()

# efficiencies (just for plotting purposes)
inFileEffD = TFile.Open("Analysis_output_LHC21dn_TOFPID.root")
hPtVsYEffD = inFileEffD.Get("hPtvsYRecSig_RecoCand")
inFileEffPi = TFile.Open("lut_eff_vs_pt.root")
gPtEffPi = inFileEffPi.Get("lutCovm.el.20kG.rmin20.geometry_v1.dat;1")
gPtEffPi.SetLineColor(kBlack)
gPtEffPi.SetLineWidth(2)

cEff = TCanvas('cEff', '', 1200, 600)
cEff.Divide(2, 1)
cEff.cd(1).SetLogz()
cEff.cd(1).SetRightMargin(0.15)
hPtVsYEffD.Draw('colz')
hPtVsYEffD.GetXaxis().SetTitleSize(0.045)
hPtVsYEffD.GetXaxis().SetLabelSize(0.045)
hPtVsYEffD.GetYaxis().SetTitleSize(0.045)
hPtVsYEffD.GetYaxis().SetLabelSize(0.045)
hPtVsYEffD.GetZaxis().SetTitleSize(0.045)
hPtVsYEffD.GetZaxis().SetLabelSize(0.045)
hPtVsYEffD.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
hPtVsYEffD.GetYaxis().SetTitle('#it{y}')
hPtVsYEffD.GetZaxis().SetTitle('D^{0} efficiency')
cEff.cd(2).DrawFrame(0., 0.5, 10., 110, ';#it{p}_{T} (GeV/#it{c});#pi efficiency (%)')
gPtEffPi.Draw('C')
gPtEffPi.GetXaxis().SetTitleSize(0.045)
gPtEffPi.GetXaxis().SetLabelSize(0.045)
gPtEffPi.GetYaxis().SetTitleSize(0.045)
gPtEffPi.GetYaxis().SetLabelSize(0.045)
gPtEffPi.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
gPtEffPi.GetYaxis().SetTitle('#pi efficiency (%)')
cEff.Modified()
cEff.Update()

predictions = {'1fm': pd.read_csv(cfg['predictions']['1fm'], names=['kstar', 'cf', 'idk'], sep=' '),
               '2fm': pd.read_csv(cfg['predictions']['2fm'], names=['kstar', 'cf', 'idk'], sep=' '),
               '3fm': pd.read_csv(cfg['predictions']['3fm'], names=['kstar', 'cf', 'idk'], sep=' '),
               '5fm': pd.read_csv(cfg['predictions']['5fm'], names=['kstar', 'cf', 'idk'], sep=' ')}

predColors = {'1fm': kAzure+4,
              '2fm': kGreen+2,
              '3fm': kOrange+7,
              '5fm': kRed+1}

gPred, sPred, hSEPred, hSEPredWithJets, hSEPredBkgD, hSEPredBkgDstar, \
    hSEPredBkgDDstar, hSEPredBkg, hCFPred = ({} for _ in range(9))
for pred in predictions:
    gPred[pred] = TGraph(1)
    gPred[pred].SetLineColor(predColors[pred])
    gPred[pred].SetFillColor(predColors[pred])
    gPred[pred].SetLineWidth(2)
    hSEPred[pred] = hSEDistr.Clone(f'hSEPred{pred}')
    hSEPred[pred].SetLineColor(predColors[pred])
    hSEPred[pred].SetMarkerColor(predColors[pred])
    for iP, (kStar, cf) in enumerate(
        zip(predictions[pred]['kstar'].to_numpy(), predictions[pred]['cf'].to_numpy())):
        gPred[pred].SetPoint(iP, kStar/1000, cf)
    sPred[pred] = TSpline3(f'sPred{pred}', gPred[pred])
    for iBin in range(1, hSEDistr.GetNbinsX()+1):
        kStarCent = hSEDistr.GetXaxis().GetBinCenter(iBin)
        hSEPred[pred].SetBinContent(
            iBin, hMEDistr.GetBinContent(iBin) * sPred[pred].Eval(kStarCent))
        hSEPred[pred].SetBinError(
            iBin, np.sqrt(hSEPred[pred].GetBinContent(iBin)))

    hSEPredBkgD[pred] = hSEPred[pred].Clone(f'hSEPredBkgD{pred}')
    hSEPredBkgDstar[pred] = hSEPred[pred].Clone(f'hSEPredBkgDstar{pred}')
    hSEPredBkgDDstar[pred] = hSEPred[pred].Clone(f'hSEPredBkgDDstar{pred}')
    hSEPredBkg[pred] = hSEPred[pred].Clone(f'hSEPredBkg{pred}')
    hSEPredWithJets[pred] = hSEPred[pred].Clone(f'hSEPredWithJets{pred}')
    hSEPredBkg[pred].SetDirectory(0)
    hSEPredWithJets[pred].SetDirectory(0)
    hSEPredBkg[pred].SetMarkerStyle(kOpenCircle)
    hSEPredWithJets[pred].SetMarkerStyle(kFullDiamond)

    for iKstarBin in range(1, hSEPred[pred].GetNbinsX()+1):
        signalSq = hSEPred[pred].GetBinContent(iKstarBin)
        signal = np.sqrt(signalSq)
        bkgD = 1./hSignalOverBkgDVsKstar.GetBinContent(iKstarBin) * signal
        bkgDstar = 1./hSignalOverBkgDstarVsKstar.GetBinContent(iKstarBin) * signal * 10
        if hSEPred[pred].GetBinContent(iBin) != np.nan:
            hSEPredBkgD[pred].SetBinContent(iKstarBin, signal * bkgD)
            hSEPredBkgDstar[pred].SetBinContent(iKstarBin, signal * bkgDstar)
            hSEPredBkgDDstar[pred].SetBinContent(iKstarBin, bkgD * bkgDstar)
            hSEPredBkg[pred].SetBinContent(
                iKstarBin, bkgD * bkgDstar + signal * bkgD + signal * bkgDstar)
        else:
            hSEPredBkg[pred].SetBinContent(1)

    hCFPred[pred] = hSEPred[pred].Clone(f'hCFPred{pred}')
    hCFPred[pred].Divide(hMEDistr)

    for iBin in range(1, hCFPred[pred].GetNbinsX()+1):
        if pred in ['1fm', '2fm']:
            S = hSEPred[pred].GetBinContent(iBin)
            B = hSEPredBkg[pred].GetBinContent(iBin)
        else:
            S = hSEPred[pred].GetBinContent(iBin) / lumiScaleFactorPP * lumiScaleFactorPbPb
            B = hSEPredBkg[pred].GetBinContent(iBin) / lumiScaleFactorPP * lumiScaleFactorPbPb

        if not isFlat:
            # we have more pairs thanks to jets, syst unc for pedestal subtraction not considered
            kStarCent = hCFPred[pred].GetBinCenter(iBin)
            S *= fPol2.Eval(kStarCent)
            hSEPredWithJets[pred].SetBinContent(iBin, S)
            if hSEPred[pred].GetBinContent(iBin) != np.nan:
                hSEPredWithJets[pred].SetBinContent(
                    iBin,
                    hSEPred[pred].GetBinContent(iBin) * fPol2.Eval(kStarCent)
                )
            else:
                hSEPredWithJets[pred].SetBinContent(1)
            B *= fPol2.Eval(kStarCent)
            hSEPredBkg[pred].SetBinContent(iBin, B)

        if hSEPred[pred].GetBinContent(iBin) != np.nan:
            hCFPred[pred].SetBinError(iBin, np.sqrt(S+2*B) / S * hCFPred[pred].GetBinContent(iBin))
        else:
            hCFPred[pred].SetBinError(iBin, 1)

legSE = TLegend(0.4, 0.2, 0.7, 0.4)
legSE.SetTextSize(0.045)
legSE.SetFillStyle(0)
legSE.SetBorderSize(0)
legSE.AddEntry(hSEPred['1fm'], 'Signal pairs', 'p')
legSE.AddEntry(hSEPredBkg['1fm'], 'Background pairs', 'p')

legModels = TLegend(0.5, 0.6, 0.85, 0.8)
legModels.SetTextSize(0.045)
legModels.SetFillStyle(0)
legModels.SetBorderSize(0)
legModels.AddEntry(hCFPred['1fm'], '1 fm (pp)', 'l')
legModels.AddEntry(hCFPred['2fm'], '2 fm', 'l')
legModels.AddEntry(hCFPred['3fm'], '3 fm', 'l')
legModels.AddEntry(hCFPred['5fm'], '5 fm (Pb#minusPb)', 'l')

cDistrPred = TCanvas('cDistrPred', '', 1250, 500)
cDistrPred.Divide(3, 1)
cDistrPred.cd(1).DrawFrame(0., 1., args.xMax, hSEPred[pred].GetMaximum()*10,
                           f';#it{{k}}* (GeV/#it{{c}});{hSEPred[pred].GetYaxis().GetTitle()}')
cDistrPred.cd(1).SetLogy()
for pred in predictions:
    if pred in ['1fm']:
        hSEPred[pred].Draw('esame')
        if not isFlat:
            hSEPredWithJets[pred].Draw('esame')
        hSEPredBkg[pred].Draw('esame')
legSE.Draw()
cDistrPred.cd(2).DrawFrame(0., 1., args.xMax, hMEDistr.GetMaximum()*10,
                           f';#it{{k}}* (GeV/#it{{c}});{hMEDistr.GetYaxis().GetTitle()}')
cDistrPred.cd(2).SetLogy()
hMEDistr.Draw('esame')
cDistrPred.cd(3).DrawFrame(0., 0., args.xMax, args.yMax,
                           f';#it{{k}}* (GeV/#it{{c}});#it{{C}}_{{{Dtitle}{Dstartitle}}}')
for pred in predictions:
    gPred[pred].Draw('C')
    if pred in ['1fm']:
        hCFPred[pred].Draw('esame')
legModels.Draw()
cDistrPred.Modified()
cDistrPred.Update()

lat = TLatex()
lat.SetNDC()
lat.SetTextSize(0.045)
lat.SetTextFont(42)
lat.SetTextColor(kBlack)

cResult = TCanvas('cResult', '', 800, 800)
hFrame = cResult.DrawFrame(0., 0.01, args.xMax, args.yMax,
                           f';#it{{k}}* (GeV/#it{{c}});#it{{C}}_{{{Dtitle}{Dstartitle}}}')
hFrame.GetYaxis().SetDecimals()
hFrame.GetXaxis().SetNdivisions(505)
for pred in predictions:
    gPred[pred].Draw('C')
    if pred in ['1fm', '5fm']:
        hCFPred[pred].Draw('esame')
lat.DrawLatex(0.18, 0.88, 'ALICE 3 upgrade projection')
lat.DrawLatex(0.18, 0.83, '|#it{y}| < 4')
lat.DrawLatex(0.4, 0.5, '#it{L}_{int}:')
lat.DrawLatex(0.4, 0.45, '   - pp = 18 fb^{#minus1}')
lat.DrawLatex(0.4, 0.4, '   - 0#minus10% Pb#minusPb = 35 nb^{#minus1}')
legModels.Draw()
cResult.Modified()
cResult.Update()

cResult.SaveAs(cfg['output']['file'])
cDistrPythia.SaveAs('DistrPythia_'+cfg['output']['file'])
cDistrPred.SaveAs('DistrPredExample_'+cfg['output']['file'])

if doAccStudy:
    hSEDistrVsYD, hSEDistrVsYDstar, cPairsVsRapidity, hKstar, hKstarRatio = {}, {}, {}, {}, {}
    cKstar = TCanvas('cKstar', '', 1200, 600)
    cKstar.Divide(2, 1)
    cKstar.cd(1).DrawFrame(0., 1.e-1, 1., 1.e5,
                           f';k* (GeV/#it{{c}}); {Dtitle}{Dstartitle} pairs (a.u.)')
    cKstar.cd(2).DrawFrame(0., 1.e-3, 1., 1.,
                           ';k* (GeV/#it{c}); ratio to ALICE3')
    legExp = TLegend(0.2, 0.7, 0.45, 0.9)
    legExp.SetTextSize(0.045)
    legExp.SetFillStyle(0)
    legExp.SetBorderSize(0)

    colors = {'ALICE3': kOrange+7, 'ALICE2': kRed+1, 'LHCb':kAzure+4}
    for det in ['ALICE3', 'ALICE2', 'LHCb']:
        if det == 'ALICE2':
            hSEDistrVsY.GetAxis(3).SetRange(1, 1)
        elif det == 'LHCb':
            hSEDistrVsY.GetAxis(3).SetRange(2, 2)
        hSEDistrVsYD[det] = hSEDistrVsY.Projection(0, 2)
        hSEDistrVsYDstar[det] = hSEDistrVsY.Projection(1, 2)
        hKstar[det] = hSEDistrVsY.Projection(2)
        hKstar[det].Rebin(5)
        hKstar[det].SetLineWidth(2)
        hKstar[det].SetMarkerStyle(kFullCircle)
        hKstar[det].SetMarkerColor(colors[det])
        hKstar[det].SetLineColor(colors[det])
        hKstarRatio[det] = hKstar[det].Clone(f'hKstarRatio{det}')
        hKstarRatio[det].Divide(hKstar[det], hKstar['ALICE3'], 1., 1., 'B')
        legExp.AddEntry(hKstar[det], det, 'lp')
        cKstar.cd(1).SetLogy()
        hKstar[det].Draw('esame')
        if det == 'LHCb':
            legExp.Draw()
        if det != 'ALICE3':
            cKstar.cd(2)
            hKstarRatio[det].Draw('esame')
        cPairsVsRapidity[det] = TCanvas(f'cPairsVsRapidity_{det}', '', 1200, 600)
        cPairsVsRapidity[det].Divide(2, 1)
        cPairsVsRapidity[det].cd(1).SetRightMargin(0.15)
        hSEDistrVsYD[det].Draw('colz')
        cPairsVsRapidity[det].cd(2).SetRightMargin(0.15)
        hSEDistrVsYDstar[det].Draw('colz')
        cPairsVsRapidity[det].SaveAs('Ydistr'+'_'+det+cfg['output']['file'].replace('ALICE3', ''))
    cKstar.Modified()
    cKstar.Update()
    cKstar.SaveAs('kStar_distr_diffExperiments'+cfg['output']['file'].replace('ALICE3', ''))

input('Press enter to exit')

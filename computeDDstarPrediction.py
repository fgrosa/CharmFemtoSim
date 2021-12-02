'''
Script for the computation of the expected D-D* correlation function
'''

import sys
import argparse
from typing import Literal
import pandas as pd
import numpy as np
import yaml
from ROOT import TFile, TGaxis, TCanvas, TF1, TGraph, TH1F, TSpline3, TLegend, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kGreen, kOrange, kRed, kBlack, kFullCircle, kOpenCircle, gStyle, kRainBow # pylint: disable=import-error,no-name-in-module

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
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlcfgFileName:
    cfg = yaml.load(ymlcfgFileName, yaml.FullLoader)

inFileNames = cfg['inputfiles']
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

nEvents = cfg['eventsperfile'] * len(inFileNames)

hSEDistrVsPt, hMEDistrVsPt, hSEDistrVsY = None, None, None
for inFileName in inFileNames:
    inFile = TFile.Open(inFileName)
    hTmpSE = inFile.Get(f'hPairSE_{Dstarspecie}_{Dspecie}')
    hTmpME = inFile.Get(f'hPairME_{Dstarspecie}_{Dspecie}')
    hTmpSE.SetDirectory(0)
    hTmpME.SetDirectory(0)
    if doAccStudy:
        hTmpVsY = inFile.Get(f'hPairVsY_{Dstarspecie}_{Dspecie}')
        hTmpVsY.SetDirectory(0)
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
hSEDistr.Rebin(25)
hMEDistr.Rebin(25)
hSEDistr.SetLineColor(kBlack)
hMEDistr.SetLineColor(kBlack)
hSEDistr.SetLineWidth(2)
hMEDistr.SetLineWidth(2)
hSEDistr.SetMarkerColor(kBlack)
hMEDistr.SetMarkerColor(kBlack)
hSEDistr.SetMarkerStyle(kFullCircle)
hMEDistr.SetMarkerStyle(kFullCircle)
hSignalOverBkgDVsKstar.Rebin(25)
hSignalOverBkgDstarVsKstar.Rebin(25)
hSignalOverBkgDVsKstar.Scale(1./25)
hSignalOverBkgDstarVsKstar.Scale(1./25)

hCFPythia = hSEDistr.Clone('hCFPythia')
hCFPythia.GetYaxis().SetDecimals()
hCFPythia.GetYaxis().SetTitle(f'#it{{C}}_{Dtitle}{Dstartitle}')
hCFPythia.Divide(hMEDistr)
fPol1 = TF1('fPol1', 'pol1', 1., 2., 1)
hCFPythia.Fit('fPol1', '0Q', '', 1., 2.)
isFlat = False
if abs(fPol1.GetParError(1) / fPol1.GetParameter(1)) < 2: #compatible with 0 within 2 sigma
    fPol0 = TF1('fPol0', 'pol0', 1., 2., 1)
    hCFPythia.Fit('fPol0', '0Q', '', 1., 2.)
    scaleFact = fPol0.GetParameter(0)
    isFlat = True # no jet background
else:
    scaleFact = (1-fPol1.GetParameter(0)) / fPol1.GetParameter(1)
hCFPythia.Scale(1./scaleFact)
hMEDistr.Scale(scaleFact)

cDistrPythia = TCanvas('cDistrPythia', '', 1800, 600)
cDistrPythia.Divide(3, 1)
cDistrPythia.cd(1).SetLogy()
hSEDistr.Draw('e')
cDistrPythia.cd(2).SetLogy()
hMEDistr.Draw('e')
cDistrPythia.cd(3)
hCFPythia.Draw('e')
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

gPred, sPred, hSEPred, hSEPredBkgD, hSEPredBkgDstar, hSEPredBkgDDstar, hSEPredBkg, hCFPred = ({} for _ in range(8))
for pred in predictions:
    gPred[pred] = TGraph(1)
    gPred[pred].SetLineColor(predColors[pred])
    gPred[pred].SetFillColor(predColors[pred])
    gPred[pred].SetLineWidth(2)
    hSEPred[pred] = hSEDistr.Clone(f'hSEPred{pred}')
    hSEPred[pred].SetLineColor(predColors[pred])
    hSEPred[pred].SetMarkerColor(predColors[pred])
    for iP, (kStar, cf) in enumerate(zip(predictions[pred]['kstar'].to_numpy(), predictions[pred]['cf'].to_numpy())):
        gPred[pred].SetPoint(iP, kStar/1000, cf)
    sPred[pred] = TSpline3(f'sPred{pred}', gPred[pred])
    for iBin in range(1, hSEDistr.GetNbinsX()+1):
        kStarCent = hSEDistr.GetXaxis().GetBinCenter(iBin)
        hSEPred[pred].SetBinContent(iBin, hMEDistr.GetBinContent(iBin) * sPred[pred].Eval(kStarCent))
        hSEPred[pred].SetBinError(iBin, np.sqrt(hSEPred[pred].GetBinContent(iBin)))

    hSEPredBkgD[pred] = hSEPred[pred].Clone(f'hSEPredBkgD{pred}')
    hSEPredBkgDstar[pred] = hSEPred[pred].Clone(f'hSEPredBkgDstar{pred}')
    hSEPredBkgDDstar[pred] = hSEPred[pred].Clone(f'hSEPredBkgDDstar{pred}')
    hSEPredBkg[pred] = hSEPred[pred].Clone(f'hSEPredBkg{pred}')
    hSEPredBkg[pred].SetMarkerStyle(kOpenCircle)

    for iKstarBin in range(1, hSEPred[pred].GetNbinsX()+1):
        signalSq = hSEPred[pred].GetBinContent(iKstarBin)
        signal = np.sqrt(signalSq)
        bkgD = 1./hSignalOverBkgDVsKstar.GetBinContent(iKstarBin) * signal
        bkgDstar = 1./hSignalOverBkgDstarVsKstar.GetBinContent(iKstarBin) * signal
        hSEPredBkgD[pred].SetBinContent(iKstarBin, signal * bkgD)
        hSEPredBkgDstar[pred].SetBinContent(iKstarBin, signal * bkgDstar)
        hSEPredBkgDDstar[pred].SetBinContent(iKstarBin, bkgD * bkgDstar)
        hSEPredBkg[pred].SetBinContent(iKstarBin, bkgD * bkgDstar + signal * bkgD + signal * bkgDstar)

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
            S *= fPol1.Eval(kStarCent)
            hCFPred[pred].SetBinError(iBin, np.sqrt(S+2*B) / S * hCFPred[pred].GetBinContent(iBin))

        hCFPred[pred].SetBinError(iBin, np.sqrt(S+2*B) / S * hCFPred[pred].GetBinContent(iBin))

legSE = TLegend(0.4, 0.2, 0.7, 0.4)
legSE.SetTextSize(0.045)
legSE.SetFillStyle(0)
legSE.SetBorderSize(0)
legSE.AddEntry(hSEPred[pred], 'Signal pairs', 'p')
legSE.AddEntry(hSEPredBkg[pred], 'Background pairs', 'p')

legModels = TLegend(0.5, 0.6, 0.85, 0.8)
legModels.SetTextSize(0.045)
legModels.SetFillStyle(0)
legModels.SetBorderSize(0)
legModels.AddEntry(hCFPred['1fm'], '1 fm (pp)', 'l')
legModels.AddEntry(hCFPred['2fm'], '2 fm', 'l')
legModels.AddEntry(hCFPred['3fm'], '3 fm', 'l')
legModels.AddEntry(hCFPred['5fm'], '5 fm (Pb#minusPb)', 'l')

cDistrPred = TCanvas('cDistrPred', '', 1800, 600)
cDistrPred.Divide(3, 1)
cDistrPred.cd(1).DrawFrame(0., 1., 2., hSEPred[pred].GetMaximum()*10,
                           f';#it{{k}}* (GeV/#it{{c}});{hSEPred[pred].GetYaxis().GetTitle()}')
cDistrPred.cd(1).SetLogy()
for pred in predictions:
    if pred in ['1fm']:
        hSEPred[pred].Draw('esame')
        hSEPredBkg[pred].Draw('esame')
legSE.Draw()
cDistrPred.cd(2).SetLogy()
hMEDistr.Draw('e')
cDistrPred.cd(3).DrawFrame(0., 0., 0.5, 5., f';#it{{k}}* (GeV/#it{{c}});#it{{C}}_{{{Dtitle}{Dstartitle}}}')
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
hFrame = cResult.DrawFrame(0., 0.01, 0.5, 8., f';#it{{k}}* (GeV/#it{{c}});#it{{C}}_{{{Dtitle}{Dstartitle}}}')
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

cResult.SaveAs(cfg['output']['file'])

if doAccStudy:
    hSEDistrVsYD = hSEDistrVsY.Project3D('xy')
    hSEDistrVsYDstar = hSEDistrVsY.Project3D('xz')
    cPairsVsRapidity = TCanvas('cPairsVsRapidity', '', 800, 400)
    cPairsVsRapidity.Divide(2, 1)
    cPairsVsRapidity.cd(1)
    hSEDistrVsYD.Draw('colz')
    cPairsVsRapidity.cd(2)
    hSEDistrVsYDstar.Draw('colz')

input('Press enter to exit')

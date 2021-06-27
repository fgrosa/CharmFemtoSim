#!/usr/bin/env python3

'''
Script to compare C(k*) of different charm-hadron species from pythia
'''

import os
import argparse
from ROOT import TFile, TCanvas, TLegend, gStyle, gROOT, kRed, kAzure, kGreen, kOrange, kBlack # pylint: disable=import-error,no-name-in-module

gStyle.SetPadTopMargin(0.035)
gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.04, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileRoot', metavar='text', default='AnalysisResults.root')
parser.add_argument('suffix', metavar='text', default='CRmode2')
parser.add_argument('--rebinPi', type=int, default=1)
parser.add_argument('--rebinK', type=int, default=1)
parser.add_argument('--rebinPr', type=int, default=1)
parser.add_argument('--log', action='store_true', default=False)
parser.add_argument('--kStarMax', type=float, default=1.)
args = parser.parse_args()

pdgCharm = [411, 421, 431, 413, 4122]
pdgLight = [211, 321, 2212]
namesLight = {211: 'pion', 321: 'kaon', 2212: 'proton'}
titlesLight = {211: '#pi^{#plus} (u#bar{d})', 321: 'K^{#plus} (u#bar{s})', 2212: 'p (uud)',
               -211: '#pi^{#minus} (#bar{u}d)', -321: 'K^{#minus} (#bar{u}s)', -2212: '#bar{p} (#bar{u}#bar{u}#bar{d})'}
titlesCharm = {411: 'D^{#plus} (c#bar{d})', 421: 'D^{0} (c#bar{u})', 431: 'D_{s}^{#plus} (c#bar{s})', 413: 'D*^{#plus} (c#bar{d})', 4122: '#Lambda_{c}^{#plus} (cud)'}
colors = {411: kGreen+2, 421: kRed+1, 431: kOrange+7, 413: kAzure+2, 4122: kBlack}
rebin = {211: args.rebinPi, 321: args.rebinK, 2212: args.rebinPr}

hPairsSE, hPairsME, hCF = {}, {}, {}

inFile = TFile.Open(args.inFileRoot)
for pdgC in pdgCharm:
    hPairsSE[pdgC], hPairsME[pdgC], hCF[pdgC] = {}, {}, {}
    for pdgL in pdgLight:
        hPairsSE[pdgC][pdgL] = {'part': inFile.Get(f'hPairSE_{pdgC}_{pdgL}'),
                                'antipart': inFile.Get(f'hPairSE_{pdgC}_{-pdgL}')}
        hPairsME[pdgC][pdgL] = {'part': inFile.Get(f'hPairME_{pdgC}_{pdgL}'),
                                'antipart': inFile.Get(f'hPairME_{pdgC}_{-pdgL}')}
        for hist in ['part', 'antipart']:
            hPairsSE[pdgC][pdgL][hist].SetDirectory(0)
            hPairsME[pdgC][pdgL][hist].SetDirectory(0)
            hPairsSE[pdgC][pdgL][hist].Sumw2()
            hPairsME[pdgC][pdgL][hist].Sumw2()
            hPairsSE[pdgC][pdgL][hist].Rebin(rebin[pdgL])
            hPairsME[pdgC][pdgL][hist].Rebin(rebin[pdgL])

        hCF[pdgC][pdgL] = {'part': hPairsSE[pdgC][pdgL]['part'].Clone(f'hCF_{pdgC}_{pdgL}'),
                           'antipart': hPairsSE[pdgC][pdgL]['antipart'].Clone(f'hCF_{pdgC}_{pdgL}')}

        for hist in ['part', 'antipart']:
            hCF[pdgC][pdgL][hist].SetDirectory(0)
            hPairsSE[pdgC][pdgL][hist].SetLineColor(colors[pdgC])
            hPairsSE[pdgC][pdgL][hist].SetMarkerColor(colors[pdgC])
            hPairsSE[pdgC][pdgL][hist].SetMarkerStyle(20)
            hPairsSE[pdgC][pdgL][hist].SetLineWidth(2)
            hPairsME[pdgC][pdgL][hist].SetLineColor(colors[pdgC])
            hPairsME[pdgC][pdgL][hist].SetMarkerColor(colors[pdgC])
            hPairsME[pdgC][pdgL][hist].SetMarkerStyle(20)
            hPairsME[pdgC][pdgL][hist].SetLineWidth(2)
            hCF[pdgC][pdgL][hist].SetLineColor(colors[pdgC])
            hCF[pdgC][pdgL][hist].SetMarkerColor(colors[pdgC])
            hCF[pdgC][pdgL][hist].SetMarkerStyle(20)
            hCF[pdgC][pdgL][hist].SetLineWidth(2)
            hCF[pdgC][pdgL][hist].Divide(hPairsSE[pdgC][pdgL][hist], hPairsME[pdgC][pdgL][hist], 1., 1./10)
            hCF[pdgC][pdgL][hist].Scale(1./hCF[pdgC][pdgL][hist].GetBinContent(hCF[pdgC][pdgL][hist].GetXaxis().FindBin(1)))
inFile.Close()

cCF, legPart, legAntiPart = {}, {}, {}
for pdgL in pdgLight:
    legPart[pdgL] = TLegend(0.5, 0.63, 0.85, 0.93)
    legPart[pdgL].SetTextSize(0.04)
    legPart[pdgL].SetBorderSize(0)
    legPart[pdgL].SetFillStyle(0)
    for pdgC in pdgCharm:
        legPart[pdgL].AddEntry(hCF[pdgC][pdgL]['part'], f'{titlesLight[pdgL]} #minus {titlesCharm[pdgC]}', 'lp')
    legAntiPart[pdgL] = TLegend(0.5, 0.63, 0.85, 0.93)
    legAntiPart[pdgL].SetTextSize(0.04)
    legAntiPart[pdgL].SetBorderSize(0)
    legAntiPart[pdgL].SetFillStyle(0)
    for pdgC in pdgCharm:
        legAntiPart[pdgL].AddEntry(hCF[pdgC][pdgL]['antipart'], f'{titlesLight[-pdgL]} #minus {titlesCharm[pdgC]}', 'lp')

    maxCF = 1000 if args.log else 2
    cCF[pdgL] = TCanvas(f'cCF_{pdgL}', '', 1200, 600)
    cCF[pdgL].Divide(2, 1)
    cCF[pdgL].cd(1).DrawFrame(0., 0.01, args.kStarMax, maxCF, ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
    if args.log:
        cCF[pdgL].cd(1).SetLogy()
    for pdgC in pdgCharm:
        hCF[pdgC][pdgL]['part'].Draw('esame')
    legPart[pdgL].Draw()
    cCF[pdgL].cd(2).DrawFrame(0., 0.01, args.kStarMax, maxCF, ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
    if args.log:
        cCF[pdgL].cd(2).SetLogy()
    for pdgC in pdgCharm:
        hCF[pdgC][pdgL]['antipart'].Draw('esame')
    legAntiPart[pdgL].Draw()

    cCF[pdgL].SaveAs(f'CF_{namesLight[pdgL]}{args.suffix}.pdf')

cCF[pdgL].Modified()
cCF[pdgL].Update()
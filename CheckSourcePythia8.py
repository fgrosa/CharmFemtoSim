'''
Script for the evaluation of the source distribution in pythia
'''

import argparse
import uproot

from ROOT import TFile, TGaxis, TCanvas, TH1F, TH2F, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kGreen, kOrange, kRed, kBlack, kGray, kMagenta, kTeal, kFullCircle, gStyle, kRainBow # pylint: disable=import-error,no-name-in-module

def ReleaseAxes(sparse, axis=-1):
    '''
    Helper function to release selections on axes of a THnSparse
    '''
    if axis >= 0:
        sparse.GetAxis(axis).SetRange(-1, -1)
    else:
        for ax in sparse.GetNdimensions():
            sparse.GetAxis(ax).SetRange(-1, -1)


gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadLeftMargin(0.13)
gStyle.SetPadRightMargin(0.035)
gStyle.SetTitleOffset(1.3, 'yz')
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(kRainBow)
TGaxis.SetMaxDigits(3)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inputFile', metavar='text', default='AnalysisResults.root')
parser.add_argument('--kStarMax', type=float, default=1.)
parser.add_argument('--mTmin', type=float, default=2.)
parser.add_argument('--mTmax', type=float, default=2.5)
parser.add_argument('--multMin', type=float, default=0.)
parser.add_argument('--multMax', type=float, default=200.)
args = parser.parse_args()

# to be set
pdgLight = [211, 321, 2212]
pdgCharm = [411]

colors = {
    211: {
        211: kBlack,
        321: kAzure+4,
        2212: kGreen+3,
        411: kRed+2,
        421: kOrange+7,
        413: kMagenta+3,
    },
    321: {
        211: kGray+2,
        321: kAzure+2,
        2212: kGreen+2,
        411: kRed,
        421: kOrange-3,
        413: kMagenta+1,
    },
    2212: {
        211: kGray,
        321: kAzure+1,
        2212: kTeal+2,
        411: kRed-7,
        421: kOrange-2,
        413: kMagenta-5,
    }
}

hSource, hSourceVsmT, hSourceVsKstar, hSourceVsMult = ({} for _ in range(4))
hSourceVsmTProfile, hSourceVsKstarProfile, hSourceVsMultProfile = ({} for _ in range(3))

tree = uproot.open(args.inputFile)['tupleRadius']
df = tree.arrays(library='pd')
df.query(f'kstar < {args.kStarMax}', inplace=True)
df.query(f'{args.multMin} < mult < {args.multMax}', inplace=True)

for pdgL in pdgLight:
    dfPart1 = df.query(f'abs(part1) == {pdgL}')
    hSource[pdgL] = {411: None, 421: None,
                     413: None, 211: None,
                     321: None, 2212: None}
    hSourceVsmT[pdgL] = {411: None, 421: None,
                         413: None, 211: None,
                         321: None, 2212: None}
    hSourceVsKstar[pdgL] = {411: None, 421: None,
                            413: None, 211: None,
                            321: None, 2212: None}
    hSourceVsMult[pdgL] = {411: None, 421: None,
                           413: None, 211: None,
                           321: None, 2212: None}
    hSourceVsmTProfile[pdgL] = {411: None, 421: None,
                                413: None, 211: None,
                                321: None, 2212: None}
    hSourceVsKstarProfile[pdgL] = {411: None, 421: None,
                                   413: None, 211: None,
                                   321: None, 2212: None}
    hSourceVsMultProfile[pdgL] = {411: None, 421: None,
                                  413: None, 211: None,
                                  321: None, 2212: None}
    for pdgC in pdgCharm:
        hSource[pdgL][pdgC] = TH1F(
            f'hSource_{pdgL}_{pdgC}',
            ';#it{r}_{0}* (fm);Normalised counts',
            100,
            0.,
            10.
        )
        hSourceVsmT[pdgL][pdgC] = TH2F(
            f'hSourceVsmT_{pdgL}_{pdgC}',
            ';#it{m}_{T} (GeV/#it{c}^{2});#it{r}_{0}* (fm)',
            100,
            1.8,
            3.8,
            100,
            0.,
            10.
        )
        hSourceVsKstar[pdgL][pdgC] = TH2F(
            f'hSourceVsKstar_{pdgL}_{pdgC}',
            ';#it{k}* (GeV/#it{c});#it{r}_{0}* (fm)',
            100,
            0.,
            1.,
            100,
            0.,
            10.
        )
        hSourceVsMult[pdgL][pdgC] = TH2F(
            f'hSourceVsMult_{pdgL}_{pdgC}',
            ';#it{N}_{ch};#it{r}_{0}* (fm)',
            100,
            -0.5,
            99.5,
            100,
            0.,
            10.
        )
        dfPart2 = dfPart1.query(f'abs(part2) == {pdgC}')
        for r, mt, kstar, mult in zip(
            dfPart2['rstar'].to_numpy(),
            dfPart2['mT'].to_numpy(),
            dfPart2['kstar'].to_numpy(),
            dfPart2['mult'].to_numpy(),
        ):
            if args.mTmin < mt < args.mTmax:
                hSource[pdgL][pdgC].Fill(r)
            hSourceVsmT[pdgL][pdgC].Fill(mt, r)
            hSourceVsKstar[pdgL][pdgC].Fill(kstar, r)
            hSourceVsMult[pdgL][pdgC].Fill(mult, r)

        hSourceVsmTProfile[pdgL][pdgC] = hSourceVsmT[pdgL][pdgC].ProfileX()
        hSourceVsKstarProfile[pdgL][pdgC] = hSourceVsKstar[pdgL][pdgC].ProfileX()
        hSourceVsMultProfile[pdgL][pdgC] = hSourceVsMult[pdgL][pdgC].ProfileX()
        hSourceVsmTProfile[pdgL][pdgC].SetLineColor(colors[pdgL][pdgC])
        hSourceVsmTProfile[pdgL][pdgC].SetMarkerColor(colors[pdgL][pdgC])
        hSourceVsmTProfile[pdgL][pdgC].SetMarkerStyle(kFullCircle)
        hSourceVsmTProfile[pdgL][pdgC].SetLineWidth(2)
        hSourceVsKstarProfile[pdgL][pdgC].SetLineColor(colors[pdgL][pdgC])
        hSourceVsKstarProfile[pdgL][pdgC].SetMarkerColor(colors[pdgL][pdgC])
        hSourceVsKstarProfile[pdgL][pdgC].SetMarkerStyle(kFullCircle)
        hSourceVsKstarProfile[pdgL][pdgC].SetLineWidth(2)
        hSourceVsMultProfile[pdgL][pdgC].SetLineColor(colors[pdgL][pdgC])
        hSourceVsMultProfile[pdgL][pdgC].SetMarkerColor(colors[pdgL][pdgC])
        hSourceVsMultProfile[pdgL][pdgC].SetMarkerStyle(kFullCircle)
        hSourceVsMultProfile[pdgL][pdgC].SetLineWidth(2)
        dfPart2 = None

    for pdgL2 in pdgLight:
        hSource[pdgL][pdgL2] = TH1F(
            f'hSource_{pdgL}_{pdgL2}',
            ';r_{0}*;Normalised counts',
            100,
            0.,
            10.
        )
        hSourceVsmT[pdgL][pdgL2] = TH2F(
            f'hSourceVsmT_{pdgL}_{pdgL2}',
            ';#it{m}_{T} (GeV/#it{c}^{2});#it{r}_{0}* (fm)',
            100,
            1.8,
            3.8,
            100,
            0.,
            10.
        )
        hSourceVsKstar[pdgL][pdgL2] = TH2F(
            f'hSourceVsKstar_{pdgL}_{pdgL2}',
            ';#it{k}* (GeV/#it{c});#it{r}_{0}* (fm)',
            100,
            0.,
            1.,
            100,
            0.,
            10.
        )
        hSourceVsMult[pdgL][pdgL2] = TH2F(
            f'hSourceVsMult_{pdgL}_{pdgL2}',
            ';#it{N}_{ch};#it{r}_{0}* (fm)',
            100,
            -0.5,
            99.5,
            100,
            0.,
            10.
        )
        dfPart2 = dfPart1.query(f'abs(part2) == {pdgL2}')
        for r, mt, kstar, mult in zip(
            dfPart2['rstar'].to_numpy(),
            dfPart2['mT'].to_numpy(),
            dfPart2['kstar'].to_numpy(),
            dfPart2['mult'].to_numpy(),
        ):
            if args.mTmin < mt < args.mTmax:
                hSource[pdgL][pdgL2].Fill(r)
            hSourceVsmT[pdgL][pdgL2].Fill(mt, r)
            hSourceVsKstar[pdgL][pdgL2].Fill(kstar, r)
            hSourceVsMult[pdgL][pdgL2].Fill(mult, r)

        hSourceVsmTProfile[pdgL][pdgL2] = hSourceVsmT[pdgL][pdgL2].ProfileX()
        hSourceVsKstarProfile[pdgL][pdgL2] = hSourceVsKstar[pdgL][pdgL2].ProfileX()
        hSourceVsMultProfile[pdgL][pdgL2] = hSourceVsMult[pdgL][pdgL2].ProfileX()
        hSourceVsmTProfile[pdgL][pdgL2].SetLineColor(colors[pdgL][pdgL2])
        hSourceVsmTProfile[pdgL][pdgL2].SetMarkerColor(colors[pdgL][pdgL2])
        hSourceVsmTProfile[pdgL][pdgL2].SetMarkerStyle(kFullCircle)
        hSourceVsmTProfile[pdgL][pdgL2].SetLineWidth(2)
        hSourceVsKstarProfile[pdgL][pdgL2].SetLineColor(colors[pdgL][pdgL2])
        hSourceVsKstarProfile[pdgL][pdgL2].SetMarkerColor(colors[pdgL][pdgL2])
        hSourceVsKstarProfile[pdgL][pdgL2].SetMarkerStyle(kFullCircle)
        hSourceVsKstarProfile[pdgL][pdgL2].SetLineWidth(2)
        hSourceVsMultProfile[pdgL][pdgL2].SetLineColor(colors[pdgL][pdgL2])
        hSourceVsMultProfile[pdgL][pdgL2].SetMarkerColor(colors[pdgL][pdgL2])
        hSourceVsMultProfile[pdgL][pdgL2].SetMarkerStyle(kFullCircle)
        hSourceVsMultProfile[pdgL][pdgL2].SetLineWidth(2)
        dfPart2 = None

    dfPart1 = None

for pdgL in pdgLight:
    for pdgC in pdgCharm:
        hSource[pdgL][pdgC].SetLineColor(colors[pdgL][pdgC])
        hSource[pdgL][pdgC].SetMarkerColor(colors[pdgL][pdgC])
        hSource[pdgL][pdgC].SetMarkerStyle(kFullCircle)
        hSource[pdgL][pdgC].SetLineWidth(2)
        if hSource[pdgL][pdgC].Integral() > 0:
            hSource[pdgL][pdgC].Scale(1./hSource[pdgL][pdgC].Integral())
    for pdgL2 in pdgLight:
        hSource[pdgL][pdgL2].SetLineColor(colors[pdgL][pdgL2])
        hSource[pdgL][pdgL2].SetMarkerColor(colors[pdgL][pdgL2])
        hSource[pdgL][pdgL2].SetMarkerStyle(kFullCircle)
        hSource[pdgL][pdgL2].SetLineWidth(2)
        if hSource[pdgL][pdgL2].Integral() > 0:
            hSource[pdgL][pdgL2].Scale(1./hSource[pdgL][pdgL2].Integral())

outFile = TFile('SourceDistr.root', 'recreate')

cSourceVsmT, cSourceVsKstar, cSourceVsMult = {}, {}, {}
for species1 in [2212, 211, 321]:
    cSourceVsmT[species1] = {}
    cSourceVsKstar[species1] = {}
    cSourceVsMult[species1] = {}
    for species2 in [2212, 211, 321, 411]:

        cSourceVsmT[species1][species2] = TCanvas(
            f'cSourceVsmT_{species1}_{species2}',
            '',
            500,
            500
        )
        hSourceVsmT[species1][species2].Draw('colz')
        cSourceVsmT[species1][species2].Modified()
        cSourceVsmT[species1][species2].Update()
        cSourceVsmT[species1][species2].SaveAs(f'SourceDistr_vs_mT_{species1}_{species2}.pdf')

        cSourceVsKstar[species1][species2] = TCanvas(
            f'cSourceVsKstar_{species1}_{species2}',
            '',
            500,
            500
        )
        hSourceVsKstar[species1][species2].Draw('colz')
        cSourceVsKstar[species1][species2].Modified()
        cSourceVsKstar[species1][species2].Update()
        cSourceVsKstar[species1][species2].SaveAs(f'SourceDistr_vs_Kstar_{species1}_{species2}.pdf')

        cSourceVsMult[species1][species2] = TCanvas(
            f'cSourceVsMult_{species1}_{species2}',
            '',
            500,
            500
        )
        hSourceVsMult[species1][species2].Draw('colz')
        cSourceVsMult[species1][species2].Modified()
        cSourceVsMult[species1][species2].Update()
        cSourceVsMult[species1][species2].SaveAs(f'SourceDistr_vs_Mult_{species1}_{species2}.pdf')

        outFile.cd()
        hSourceVsmT[species1][species2].Write()
        hSourceVsKstar[species1][species2].Write()
        hSourceVsMult[species1][species2].Write()

leg = TLegend(0.5, 0.6, 0.8, 0.9)
leg.SetTextSize(0.04)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(hSource[2212][2212], 'pp', 'p')
leg.AddEntry(hSource[321][321], 'KK', 'p')
leg.AddEntry(hSource[211][211], '#pi#pi', 'p')
leg.AddEntry(hSource[2212][411], 'pD', 'p')
leg.AddEntry(hSource[211][411], '#piD', 'p')
leg.AddEntry(hSource[321][411], 'KD', 'p')

legWmean = TLegend(0.5, 0.6, 0.8, 0.9)
legWmean.SetTextSize(0.04)
legWmean.SetFillStyle(0)
legWmean.SetBorderSize(0)
legWmean.AddEntry(hSource[2212][2212], f'pp, #LT #it{{r}}_{{0}}* #GT = {hSource[2212][2212].GetMean():0.2f}#pm{hSource[2212][2212].GetMeanError():0.2f}', 'p')
legWmean.AddEntry(hSource[321][321], f'KK, #LT #it{{r}}_{{0}}* #GT = {hSource[321][321].GetMean():0.2f}#pm{hSource[321][321].GetMeanError():0.2f}', 'p')
legWmean.AddEntry(hSource[211][211], f'#pi#pi, #LT #it{{r}}_{{0}}* #GT = {hSource[211][211].GetMean():0.2f}#pm{hSource[211][211].GetMeanError():0.2f}', 'p')
legWmean.AddEntry(hSource[2212][411], f'pD, #LT #it{{r}}_{{0}}* #GT = {hSource[2212][411].GetMean():0.2f}#pm{hSource[2212][411].GetMeanError():0.2f}', 'p')
legWmean.AddEntry(hSource[211][411], f'#piD, #LT #it{{r}}_{{0}}* #GT = {hSource[211][411].GetMean():0.2f}#pm{hSource[211][411].GetMeanError():0.2f}', 'p')
legWmean.AddEntry(hSource[321][411], f'KD, #LT #it{{r}}_{{0}}* #GT = {hSource[321][411].GetMean():0.2f}#pm{hSource[321][411].GetMeanError():0.2f}', 'p')

cSourceVsmTProfile = TCanvas('cSourceVsmTProfile', '', 500, 500)
hFrame = cSourceVsmTProfile.DrawFrame(1.8, 0.1, 3.8, 7., ';#it{m}_{T} (GeV/#it{c}^{2});#LT #it{r}_{0}* #GT (fm)')
hFrame.GetYaxis().SetDecimals()
hSourceVsmTProfile[2212][2212].Draw('esame')
hSourceVsmTProfile[321][321].Draw('esame')
hSourceVsmTProfile[211][211].Draw('esame')
hSourceVsmTProfile[2212][411].Draw('esame')
hSourceVsmTProfile[211][411].Draw('esame')
hSourceVsmTProfile[321][411].Draw('esame')
leg.Draw()
cSourceVsmTProfile.Modified()
cSourceVsmTProfile.Update()
outFile.cd()
for species1 in [2212, 211, 321]:
    for species2 in [2212, 411]:
        hSourceVsmTProfile[species1][species2].Write()
cSourceVsmTProfile.Write()
cSourceVsmTProfile.SaveAs('cSourceVsmTProfile.pdf')

cSourceVsKstarProfile = TCanvas('cSourceVsKstarProfile', '', 500, 500)
hFrame = cSourceVsKstarProfile.DrawFrame(0., 0.1, 1., 7., ';#it{k}* (GeV/#it{c});#LT #it{r}_{0}* #GT (fm)')
hFrame.GetYaxis().SetDecimals()
hSourceVsKstarProfile[2212][2212].Draw('esame')
hSourceVsKstarProfile[321][321].Draw('esame')
hSourceVsKstarProfile[211][211].Draw('esame')
hSourceVsKstarProfile[2212][411].Draw('esame')
hSourceVsKstarProfile[211][411].Draw('esame')
hSourceVsKstarProfile[321][411].Draw('esame')
leg.Draw()
cSourceVsKstarProfile.Modified()
cSourceVsKstarProfile.Update()
outFile.cd()
for species1 in [2212, 211, 321]:
    for species2 in [2212, 411]:
        hSourceVsKstarProfile[species1][species2].Write()
cSourceVsKstarProfile.Write()
cSourceVsKstarProfile.SaveAs('cSourceVsKstarProfile.pdf')

cSourceVsMultProfile = TCanvas('cSourceVsMultProfile', '', 500, 500)
hFrame = cSourceVsMultProfile.DrawFrame(-0.5, 0.1, 49.5, 7., ';#it{N}_{ch};#LT #it{r}_{0}* #GT (fm)')
hFrame.GetYaxis().SetDecimals()
hSourceVsMultProfile[2212][2212].Draw('esame')
hSourceVsMultProfile[321][321].Draw('esame')
hSourceVsMultProfile[211][211].Draw('esame')
hSourceVsMultProfile[2212][411].Draw('esame')
hSourceVsMultProfile[211][411].Draw('esame')
hSourceVsMultProfile[321][411].Draw('esame')
leg.Draw()
cSourceVsMultProfile.Modified()
cSourceVsMultProfile.Update()
outFile.cd()
for species1 in [2212, 211, 321]:
    for species2 in [2212, 411]:
        hSourceVsMultProfile[species1][species2].Write()
cSourceVsMultProfile.Write()
cSourceVsMultProfile.SaveAs('cSourceVsMultProfile.pdf')

cSourceDistr = TCanvas('cSourceDistr', '', 500, 500)
cSourceDistr.DrawFrame(0., 1.e-4, 6., 1., ';#it{r}_{0}* (fm);Normalised counts')
cSourceDistr.SetLogy()
hSource[2212][2212].Draw('esame')
hSource[321][321].Draw('esame')
hSource[211][211].Draw('esame')
hSource[2212][411].Draw('esame')
hSource[211][411].Draw('esame')
hSource[321][411].Draw('esame')
legWmean.Draw()
cSourceDistr.Modified()
cSourceDistr.Update()
outFile.cd()
for species1 in [2212, 211, 321]:
    for species2 in [2212, 411]:
        hSource[species1][species2].SetNameTitle(
            f'hSource_{species1}_{species2}',
            f';r_{{0}}*({species1}, {species2}) (fm);Normalised counts'
        )
        hSource[species1][species2].Write()
cSourceDistr.Write()

cSourceDistr.SaveAs('SourceDistr.pdf')

input()

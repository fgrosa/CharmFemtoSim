#if !defined(__CINT__) || defined(__MAKECINT__)

#include <array>
#include <string>
#include <vector>
#include <map>
#include <deque>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TSpline.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"

#endif

using namespace Pythia8;

namespace
{
    enum tunes
    {
        kMonash = 0,
        kCRMode0,
        kCRMode2,
        kCRMode3
    };

    enum processes
    {
        kSoftQCD = 0,
        kHardQCD
    };

    std::array<int, 1> DmesonPDG{421}; // D0
    std::array<int, 1> DstarPDG{413}; // D*+
}

//__________________________________________________________________________________________________
void SimulateDDstarCorrelation(int nEvents=1000000, int tune=kCRMode2, int process=kSoftQCD, float energy=13000, int seed=42, std::string outFileNameRoot="AnalysisResults.root");
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2);
template<typename T>
bool CheckDauAcc(T &ptDauArray, T &etaDauArray);
template<typename T, typename T2>
bool IsFromBeauty(T &mothers, T2 &pythia);

//__________________________________________________________________________________________________
void SimulateDDstarCorrelation(int nEvents, int tune, int process, float energy, int seed, std::string outFileNameRoot)
{
    //__________________________________________________________
    // create and configure pythia generator

    Pythia pythia;
    if(process == kSoftQCD)
    {
        pythia.readString("SoftQCD:all = on");
    }
    else if(process == kHardQCD)
    {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // set tune
    if(tune == kMonash)
    {
        pythia.readString(Form("Tune:pp = 14"));
    }
    else if(tune == kCRMode0)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode2)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode3)
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }

    // keep only interesting decays, to be reweighted a posteriori
    pythia.readString("421:onMode = off");
    pythia.readString("411:onMode = off");
    pythia.readString("413:onMode = off");
    pythia.readString("423:onMode = off");
    pythia.readString("421:onIfMatch = 211 321");
    pythia.readString("411:onIfMatch = 211 211 321");
    pythia.readString("413:onIfMatch = 211 421");
    pythia.readString("423:onIfMatch = 11 421");

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    // load efficiencies and mass distributions
    TFile* inFileEffD = TFile::Open("output_PbPbminbias.root"); //PbPb file
    TH1F* hPtEffD[4];
    for(auto iRap=0; iRap<4; iRap++)
    {
        hPtEffD[iRap] = (TH1F*)inFileEffD->Get(Form("Efficiency_TOFplusRICHPID_rap%d", iRap));
    }
    TFile* inFileEffPi = TFile::Open("lut_eff_vs_pt.root");
    TGraph* gPtEffPi = (TGraph*)inFileEffPi->Get("lutCovm.el.20kG.rmin20.geometry_v1.dat;1");
    TSpline3* sPtEffPi = new TSpline3("sPtEffPi", gPtEffPi);

    //__________________________________________________________
    // define outputs
    std::map<int, std::map<int, std::map<std::string, TH3F*>>> hPairSE, hPairME; // all combinations of D, D*, particle, antiparticle
    for(auto &pdgDstar: DstarPDG)
    {
        for(auto &pdgDmeson: DmesonPDG)
        {
            hPairSE[pdgDstar][pdgDmeson]["part"] = new TH3F(Form("hPairSE_%d_%d", pdgDstar, pdgDmeson), ";#it{p}_{T}^{D} (GeV/#it{c});;#it{p}_{T}^{D^{*}} (GeV/#it{c});#it{k}^{*} (GeV/#it{c})pairs", 100, 0., 10., 100, 0., 10., 2000, 0., 2.);
            hPairSE[pdgDstar][pdgDmeson]["antipart"] = new TH3F(Form("hPairSE_%d_%d", pdgDstar, -pdgDmeson), ";#it{p}_{T}^{D} (GeV/#it{c});;#it{p}_{T}^{D^{*}} (GeV/#it{c});#it{k}^{*} (GeV/#it{c});pairs", 100, 0., 10., 100, 0., 10., 2000, 0., 2.);
            hPairME[pdgDstar][pdgDmeson]["part"] = new TH3F(Form("hPairME_%d_%d", pdgDstar, pdgDmeson), ";#it{p}_{T}^{D} (GeV/#it{c});;#it{p}_{T}^{D^{*}} (GeV/#it{c});#it{k}^{*} (GeV/#it{c});pairs", 100, 0., 10., 100, 0., 10., 2000, 0., 2.);
            hPairME[pdgDstar][pdgDmeson]["antipart"] = new TH3F(Form("hPairME_%d_%d", pdgDstar, -pdgDmeson), ";#it{p}_{T}^{D} (GeV/#it{c});;#it{p}_{T}^{D^{*}} (GeV/#it{c});#it{k}^{*} (GeV/#it{c});pairs", 100, 0., 10., 100, 0., 10., 2000, 0., 2.);
        }
    }

    //__________________________________________________________
    // perform the simulation
    std::vector<ROOT::Math::PxPyPzMVector> partDmeson{};
    std::vector<ROOT::Math::PxPyPzMVector> partDstar{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferDmeson{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferDstar{};
    std::vector<int> pdgDmeson{};
    std::vector<int> pdgDstar{};
    std::vector<float> effDmeson{};
    std::vector<float> effDstar{};
    std::deque<std::vector<int>> pdgBufferDmeson{};
    std::deque<std::vector<int>> pdgBufferDstar{};
    std::deque<std::vector<float>> effBufferDmeson{};
    std::deque<std::vector<float>> effBufferDstar{};

    for (int iEvent=0; iEvent<nEvents; iEvent++)
    {
        pythia.next();
        for(int iPart=3; iPart<pythia.event.size(); iPart++)
        {
            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);

            if (find(DmesonPDG.begin(), DmesonPDG.end(), absPdg) == DmesonPDG.end() && find(DstarPDG.begin(), DstarPDG.end(), absPdg) == DstarPDG.end())
                continue;

            std::vector<int> mothers = pythia.event[iPart].motherList();
            bool isFromDstar = false;
            if(absPdg == 421) // check if D0 is D* daughter, reject
            {
                for(auto &mom: mothers)
                {
                    if(std::abs(pythia.event[mom].id()) == 413) // like setting a veto
                        isFromDstar = true;
                }
            }
            if(isFromDstar)
                continue;

            auto dauList = pythia.event[iPart].daughterList();
            std::vector<double> ptDau{}, etaDau{};
            for(auto &dau: dauList)
            {
                auto pdgDau = std::abs(pythia.event[iPart].id());
                if((absPdg != 423 && (pdgDau == 211 || pdgDau == 321)) ||
                   (absPdg == 423 && (pdgDau == 211 || pdgDau == 321 || pdgDau == 22)))
                {
                    ptDau.push_back(std::sqrt(pythia.event[dau].px()*pythia.event[dau].px() + pythia.event[dau].py()*pythia.event[dau].py() + pythia.event[dau].pz()*pythia.event[dau].pz()));
                    etaDau.push_back(pythia.event[dau].eta());
                }
            }

            if(std::find(DstarPDG.begin(), DstarPDG.end(), absPdg) != DstarPDG.end())
            {
                bool isFromB = IsFromBeauty(mothers, pythia);
                if(!isFromB)
                {
                    if(!CheckDauAcc(ptDau, etaDau))
                        continue;

                    auto ptD = std::sqrt(pythia.event[dauList[0]].px()*pythia.event[dauList[0]].px() + pythia.event[dauList[0]].py()*pythia.event[dauList[0]].py());
                    auto ptPi = std::sqrt(pythia.event[dauList[1]].px()*pythia.event[dauList[1]].px() + pythia.event[dauList[1]].py()*pythia.event[dauList[1]].py());

                    auto yPart = abs(pythia.event[dauList[0]].y());
                    int iY = -1;
                    if(yPart>4)
                        continue;
                    else if(yPart<1)
                        iY = 0;
                    else if(yPart<2)
                        iY = 1;
                    else if(yPart<3)
                        iY = 2;
                    else
                        iY = 3;
                    auto binPt = hPtEffD[iY]->GetXaxis()->FindBin(ptD);
                    auto effD = hPtEffD[iY]->GetBinContent(binPt);

                    auto effPi = sPtEffPi->Eval(ptPi) / 100;
                    effDstar.push_back(effD*effPi);

                    ROOT::Math::PxPyPzMVector part(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                    partDstar.push_back(part);
                    pdgDstar.push_back(pdg);
                }
            }
            else if(std::find(DmesonPDG.begin(), DmesonPDG.end(), absPdg) != DmesonPDG.end())
            {
                bool isFromB = IsFromBeauty(mothers, pythia);
                if(!isFromB)
                {
                    if(!CheckDauAcc(ptDau, etaDau))
                        continue;

                    auto ptPart = std::sqrt(pythia.event[iPart].px()*pythia.event[iPart].px() + pythia.event[iPart].py()*pythia.event[iPart].py());
                    auto yPart = abs(pythia.event[iPart].y());
                    int iY = -1;
                    if(yPart>4)
                        continue;
                    else if(yPart<1)
                        iY = 0;
                    else if(yPart<2)
                        iY = 1;
                    else if(yPart<3)
                        iY = 2;
                    else
                        iY = 3;
                    auto binPt = hPtEffD[iY]->GetXaxis()->FindBin(ptPart);
                    effDmeson.push_back(hPtEffD[iY]->GetBinContent(binPt));

                    ROOT::Math::PxPyPzMVector part(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                    partDmeson.push_back(part);
                    pdgDmeson.push_back(pdg);
                }
            }
        }

        if(partDmeson.size() > 0)
        {
            partBufferDmeson.push_back(partDmeson);
            pdgBufferDmeson.push_back(pdgDmeson);
            effBufferDmeson.push_back(effDmeson);
        }
        if(partDstar.size() > 0)
        {
            partBufferDstar.push_back(partDstar);
            pdgBufferDstar.push_back(pdgDstar);
            effBufferDstar.push_back(effDstar);
        }
        if (partBufferDmeson.size() > 10) { //buffer full, let's kill the first entry
            partBufferDmeson.pop_front();
            pdgBufferDmeson.pop_front();
            effBufferDmeson.pop_front();
        }
        if (partBufferDmeson.size() > 10) { //buffer full, let's kill the first entry
            partBufferDstar.pop_front();
            pdgBufferDstar.pop_front();
            effBufferDstar.pop_front();
        }

        // same event
        for(size_t iDstar=0; iDstar<partDstar.size(); iDstar++)
        {
            for(size_t iDmeson=0; iDmeson<partDmeson.size(); iDmeson++)
            {
                double kStar = ComputeKstar(partDstar[iDstar], partDmeson[iDmeson]);
                double pTDstar = partDstar[iDstar].pt();
                double pTDmeson = partDmeson[iDmeson].pt();
                if(pdgDstar[iDstar] * pdgDmeson[iDmeson] > 0)
                    hPairSE[std::abs(pdgDstar[iDstar])][std::abs(pdgDmeson[iDmeson])]["part"]->Fill(pTDmeson, pTDstar, kStar, effDstar[iDstar]*effDmeson[iDmeson]);
                else
                    hPairSE[std::abs(pdgDstar[iDstar])][std::abs(pdgDmeson[iDmeson])]["antipart"]->Fill(pTDmeson, pTDstar, kStar, effDstar[iDstar]*effDmeson[iDmeson]);
            }
        }

        // mixed event
        if(partBufferDmeson.size() < 2 || partBufferDstar.size() < 2) // to avoid repetitions
            continue;

        for(size_t iDstar=0; iDstar<partBufferDstar[partBufferDstar.size()-1].size(); iDstar++) // last only
        {
            for(size_t iME=0; iME<partBufferDmeson.size()-1; iME++) // from 0 to last-1
            {
                for(size_t iDmeson=0; iDmeson<partBufferDmeson[iME].size(); iDmeson++)
                {
                    double kStar = ComputeKstar(partBufferDstar[partBufferDstar.size()-1][iDstar], partBufferDmeson[iME][iDmeson]);
                    double pTDstar = partDstar[iDstar].pt();
                    double pTDmeson = partDmeson[iDmeson].pt();
                    if(pdgBufferDstar[partBufferDstar.size()-1][iDstar] * pdgBufferDmeson[iME][iDmeson] > 0)
                        hPairME[std::abs(pdgBufferDstar[partBufferDstar.size()-1][iDstar])][std::abs(pdgBufferDmeson[iME][iDmeson])]["part"]->Fill(pTDmeson, pTDstar, kStar, effBufferDstar[partBufferDstar.size()-1][iDstar]*effBufferDmeson[iME][iDmeson]);
                    else
                        hPairME[std::abs(pdgBufferDstar[partBufferDstar.size()-1][iDstar])][std::abs(pdgBufferDmeson[iME][iDmeson])]["antipart"]->Fill(pTDmeson, pTDstar, kStar, effBufferDstar[partBufferDstar.size()-1][iDstar]*effBufferDmeson[iME][iDmeson]);
                }
            }
        }

        partDstar.clear();
        pdgDstar.clear();
        effDstar.clear();
        partDmeson.clear();
        pdgDmeson.clear();
        effDmeson.clear();
    }

    // save root output file
    TFile outFile(outFileNameRoot.data(), "recreate");
    for(auto &pdgDstar: DstarPDG)
    {
        for(auto &pdgDmeson: DmesonPDG)
        {
            hPairSE[pdgDstar][pdgDmeson]["part"]->Write();
            hPairSE[pdgDstar][pdgDmeson]["antipart"]->Write();
            hPairME[pdgDstar][pdgDmeson]["part"]->Write();
            hPairME[pdgDstar][pdgDmeson]["antipart"]->Write();
        }
    }
    outFile.Close();
}

//__________________________________________________________________________________________________
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2)
{
    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);

    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;
    float kStar = 0.5 * trackRelK.P();
    return kStar;
}

//__________________________________________________________________________________________________
template<typename T>
bool CheckDauAcc(T &ptDauArray, T &etaDauArray)
{
    for(size_t iDau=0; iDau<ptDauArray.size(); iDau++)
    {
        if(ptDauArray[iDau] < 0.005 || std::abs(etaDauArray[iDau]) > 4.) // pT=50 MeV and eta=4
            return false;
    }

    return true;
}

//__________________________________________________________________________________________________
template<typename T, typename T2>
bool IsFromBeauty(T &mothers, T2 &pythia)
{
    for(auto &mom: mothers)
    {
        int absPdgMom = std::abs(pythia.event[mom].id());
        if(absPdgMom == 5 || absPdgMom/100 == 5 || absPdgMom/1000 == 5 ||
           (absPdgMom-10000)/100 == 5 || (absPdgMom-20000)/100 == 5 || (absPdgMom-30000)/100 == 5 ||
           (absPdgMom-100000)/100 == 5 || (absPdgMom-200000)/100 == 5 || (absPdgMom-300000)/100 == 5)
        {  // remove beauty feed-down
            return true;
        }
    }
    return false;
}
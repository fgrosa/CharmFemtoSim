#define __ALI__ //__O2__ or __ALI__

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
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#ifdef __O2__
#include "Pythia8/Pythia.h"
#endif

#ifdef __ALI__
#include <TClonesArray.h>
#include <TParticle.h>
#include "AliPythia8.h"
#endif

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

    std::array<int, 3> lightPDG{211, 321, 2212}; // pions, kaons, protons
    std::array<int, 5> charmPDG{411, 421, 431, 413, 4122}; // D+, D0, Ds, D*, Lc
}

//__________________________________________________________________________________________________
void SimulateCharmLightCorrelation(int nEvents=1000000, int tune=kCRMode2, int process=kSoftQCD, float energy=13000, int seed=42, std::string outFileNameRoot="AnalysisResults.root");
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2);

//__________________________________________________________________________________________________
void SimulateCharmLightCorrelation(int nEvents, int tune, int process, float energy, int seed, std::string outFileNameRoot)
{
    //__________________________________________________________
    // create and configure pythia generator

#ifdef __O2__

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

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();
#endif

#ifdef __ALI__
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    AliPythia8 pythia;
    if(process == kSoftQCD)
    {
        pythia.ReadString("SoftQCD:all = on");
    }
    else if(process == kHardQCD)
    {
        pythia.ReadString("HardQCD:hardccbar = on");
        pythia.ReadString("HardQCD:hardbbbar = on");
    }

    // set tune
    if(tune == kMonash)
    {
        pythia.ReadString(Form("Tune:pp = 14"));
    }
    else if(tune == kCRMode0)
    {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 2.9");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.43");
        pythia.ReadString("ColourReconnection:timeDilationMode = 0");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode2)
    {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 0.3");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.20");
        pythia.ReadString("ColourReconnection:timeDilationMode = 2");
        pythia.ReadString("ColourReconnection:timeDilationPar = 0.18");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    }
    else if(tune == kCRMode3)
    {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 0.3");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.15");
        pythia.ReadString("ColourReconnection:timeDilationMode = 3");
        pythia.ReadString("ColourReconnection:timeDilationPar = 0.073");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    }

    // init
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString(Form("Random:seed %d", seed));
    pythia.Initialize(2212, 2212, energy);
#endif

    //__________________________________________________________
    // define outputs
    std::map<int, std::map<int, std::map<std::string, TH1F*>>> hPairSE, hPairME; // all combinations of charm, light hadrons and particle/antiparticle
    for(auto &pdgC: charmPDG)
    {
        for(auto &pdgL: lightPDG)
        {
            hPairSE[pdgC][pdgL]["part"] = new TH1F(Form("hPairSE_%d_%d", pdgC, pdgL), ";#it{k}^{*} (GeV/#it{c});pairs", 2000, 0., 2.);
            hPairSE[pdgC][pdgL]["antipart"] = new TH1F(Form("hPairSE_%d_%d", pdgC, -pdgL), ";#it{k}^{*} (GeV/#it{c});pairs", 2000, 0., 2.);
            hPairME[pdgC][pdgL]["part"] = new TH1F(Form("hPairME_%d_%d", pdgC, pdgL), ";#it{k}^{*} (GeV/#it{c});pairs", 2000, 0., 2.);
            hPairME[pdgC][pdgL]["antipart"] = new TH1F(Form("hPairME_%d_%d", pdgC, -pdgL), ";#it{k}^{*} (GeV/#it{c});pairs", 2000, 0., 2.);
        }
    }

    //__________________________________________________________
    // perform the simulation
    std::vector<ROOT::Math::PxPyPzMVector> partLF{};
    std::vector<ROOT::Math::PxPyPzMVector> partCharm{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferLF{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferCharm{};
    std::vector<int> idxCharm{};
    std::vector<std::vector<int>> motherLF{};
    std::vector<int> idxLF{};
    std::vector<int> pdgLF{};
    std::vector<int> pdgCharm{};
    std::deque<std::vector<int>> pdgBufferLF{};
    std::deque<std::vector<int>> pdgBufferCharm{};

#ifdef __ALI__
TClonesArray* particles = new TClonesArray("TParticle", 1000);
#endif

    for (int iEvent=0; iEvent<nEvents; iEvent++)
    {         
       
#ifdef __O2__
        pythia.next();
        for(int iPart=3; iPart<pythia.event.size(); iPart++)
        {
            if(std::abs(pythia.event[iPart].y()) > 2) // keep only midrapidity
                continue;

            double prodR = std::sqrt(pythia.event[iPart].xProd()*pythia.event[iPart].xProd() + 
                                     pythia.event[iPart].yProd()*pythia.event[iPart].yProd() + 
                                     pythia.event[iPart].zProd()*pythia.event[iPart].zProd());

            if(prodR > 1.) // to remove products of strange decays 
                continue;

            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);
            if(std::find(charmPDG.begin(), charmPDG.end(), absPdg) != charmPDG.end())
            {
                bool isFromB = false;
                std::vector<int> mothers = pythia.event[iPart].motherList();
                for(auto &mom: mothers)
                {
                    int absPdgMom = std::abs(pythia.event[mom].id());
                    if(absPdgMom == 5 || absPdgMom/100 == 5 || absPdgMom/1000 == 5 ||
                       (absPdgMom-10000)/100 == 5 || (absPdgMom-20000)/100 == 5 || (absPdgMom-30000)/100 == 5 ||
                       (absPdgMom-100000)/100 == 5 || (absPdgMom-200000)/100 == 5 || (absPdgMom-300000)/100 == 5)
                    {  // remove beauty feed-down
                        isFromB = true;
                        break;
                    }
                }
                if(!isFromB)
                {
                    ROOT::Math::PxPyPzMVector part(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                    partCharm.push_back(part);
                    pdgCharm.push_back(pdg);
                    idxCharm.push_back(iPart);
                }
            }
            else if(std::find(lightPDG.begin(), lightPDG.end(), absPdg) != lightPDG.end())
            {
                std::vector<int> mothers = pythia.event[iPart].motherList();
                motherLF.push_back(mothers);
                ROOT::Math::PxPyPzMVector part(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                partLF.push_back(part);
                pdgLF.push_back(pdg);
            }
        }
#endif

#ifdef __ALI__
        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");
        for(int iPart=2; iPart<particles->GetEntriesFast(); iPart++)
        {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
            if(std::abs(particle->Y()) > 2) // keep only midrapidity
                continue;

            double prodR = std::sqrt(particle->Vx()*particle->Vx() + 
                                     particle->Vy()*particle->Vy() + 
                                     particle->Vz()*particle->Vz());

            if(prodR > 1.) // to remove products of strange decays 
                continue;

            int pdg = particle->GetPdgCode();
            int absPdg = std::abs(pdg);
            if(std::find(charmPDG.begin(), charmPDG.end(), absPdg) != charmPDG.end())
            {
                bool isFromB = false;
                int motherIdx = particle->GetFirstMother();
                while(motherIdx>1) // 0 and 1 protons
                {
                    int absPdgMom = std::abs(particle->GetPdgCode());
                    if(absPdgMom == 5 || absPdgMom/100 == 5 || absPdgMom/1000 == 5 ||
                       (absPdgMom-10000)/100 == 5 || (absPdgMom-20000)/100 == 5 || (absPdgMom-30000)/100 == 5 ||
                       (absPdgMom-100000)/100 == 5 || (absPdgMom-200000)/100 == 5 || (absPdgMom-300000)/100 == 5)
                    {  // remove beauty feed-down
                        isFromB = true;
                        break;
                    }
                    TParticle* mom = dynamic_cast<TParticle*>(particles->At(motherIdx));
                    motherIdx = mom->GetFirstMother();
                }
                if(!isFromB)
                {
                    ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                    partCharm.push_back(part);
                    pdgCharm.push_back(pdg);
                    idxCharm.push_back(iPart);
                }
            }
            else if(std::find(lightPDG.begin(), lightPDG.end(), absPdg) != lightPDG.end())
            {
                std::vector<int> mothers{};
                int motherIdx = particle->GetFirstMother();
                while(motherIdx>1) // 0 and 1 protons
                {
                    mothers.push_back(motherIdx);
                    TParticle* mom = dynamic_cast<TParticle*>(particles->At(motherIdx));
                    motherIdx = mom->GetFirstMother();
                }
                motherLF.push_back(mothers);
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
                partLF.push_back(part);
                pdgLF.push_back(pdg);
            }
        }
#endif

        partBufferLF.push_back(partLF);
        partBufferCharm.push_back(partCharm);
        pdgBufferLF.push_back(pdgLF);
        pdgBufferCharm.push_back(pdgCharm);
        if (partBufferLF.size() > 10) { //buffer full, let's kill the first entry
            partBufferLF.pop_front();
            partBufferCharm.pop_front();
            pdgBufferLF.pop_front();
            pdgBufferCharm.pop_front();
        } 

        // same event
        for(size_t iCharm=0; iCharm<partCharm.size(); iCharm++)
        {            
            for(size_t iLF=0; iLF<partLF.size(); iLF++)
            {
                if(std::find(motherLF[iLF].begin(), motherLF[iLF].end(), idxCharm[iCharm]) != motherLF[iLF].end())
                    continue;

                double kStar = ComputeKstar(partCharm[iCharm], partLF[iLF]);
                if(pdgCharm[iCharm] * pdgLF[iLF] > 0)
                    hPairSE[std::abs(pdgCharm[iCharm])][std::abs(pdgLF[iLF])]["part"]->Fill(kStar);
                else
                    hPairSE[std::abs(pdgCharm[iCharm])][std::abs(pdgLF[iLF])]["antipart"]->Fill(kStar);
            }
        }

        // mixed event
        if(partBufferLF.size() < 2) // to avoid repetitions
            continue;

        for(size_t iCharm=0; iCharm<partBufferCharm[partBufferCharm.size()-1].size(); iCharm++) // last only
        {            
            for(size_t iME=0; iME<partBufferLF.size()-1; iME++) // from 0 to last-1
            {
                for(size_t iLF=0; iLF<partBufferLF[iME].size(); iLF++)
                {
                    double kStar = ComputeKstar(partBufferCharm[partBufferCharm.size()-1][iCharm], partBufferLF[iME][iLF]);
                    if(pdgBufferCharm[partBufferCharm.size()-1][iCharm] * pdgBufferLF[iME][iLF] > 0)
                        hPairME[std::abs(pdgBufferCharm[partBufferCharm.size()-1][iCharm])][std::abs(pdgBufferLF[iME][iLF])]["part"]->Fill(kStar);
                    else
                        hPairME[std::abs(pdgBufferCharm[partBufferCharm.size()-1][iCharm])][std::abs(pdgBufferLF[iME][iLF])]["antipart"]->Fill(kStar);
                }
            }
        }

        partCharm.clear();
        pdgCharm.clear();
        idxCharm.clear();
        partLF.clear();
        pdgLF.clear();
        motherLF.clear();
    }

    // save root output file
    TFile outFile(outFileNameRoot.data(), "recreate");
    for(auto &pdgC: charmPDG)
    {
        for(auto &pdgL: lightPDG)
        {
            hPairSE[pdgC][pdgL]["part"]->Write();
            hPairSE[pdgC][pdgL]["antipart"]->Write();
            hPairME[pdgC][pdgL]["part"]->Write();
            hPairME[pdgC][pdgL]["antipart"]->Write();
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
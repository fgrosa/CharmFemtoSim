#if !defined(__CINT__) || defined(__MAKECINT__)

#include <array>
#include <string>
#include <vector>
#include <map>
#include <deque>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TNtuple.h>
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

//__________________________________________________________________________________________________
void SimulateSourceCharmHadrons(int nEvents=1e6, int seed=42, std::string outFileName="AnalsysiResults.root", std::string tune="CRmode2");
void ComputeRandKstar(ROOT::Math::PxPyPzMVector mom1, ROOT::Math::PxPyPzMVector mom2, ROOT::Math::XYZTVector vtx1, ROOT::Math::XYZTVector vtx2, float &kStar, float &rStar, float &mT);
bool IsStable(int absPdg);
bool IsPrimary(int absPdg);

//__________________________________________________________________________________________________
void SimulateSourceCharmHadrons(int nEvents, int seed, std::string outFileName, std::string tune) {

    //__________________________________________________________
    // create and configure pythia generator
    Pythia pythia;
    pythia.readString("SoftQCD:all = on");
    pythia.readString("Tune:pp = 14");
    if(tune == "CRmode2") {
        // set tune (CR mode2)
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
    // needed for production vertex
    pythia.readString("Fragmentation:setVertices = on");
    pythia.readString("PartonVertex:setVertex = on");

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", 13000.); //GeV
    pythia.init();

    std::vector<int> pdgCharm{};
    std::vector<int> pdgLight{};
    std::vector<ROOT::Math::PxPyPzMVector> momCharm{};
    std::vector<ROOT::Math::PxPyPzMVector> momLight{};
    std::vector<ROOT::Math::XYZTVector> prodvtxCharm{};
    std::vector<ROOT::Math::XYZTVector> prodvtxLight{};
    std::vector<int> pdgMotherCharm{};
    std::vector<int> pdgMotherLight{};
    std::vector<int> idxMotherCharm{};
    std::vector<int> idxMotherLight{};

    TNtuple* tupleRadius = new TNtuple("tupleRadius", "tupleRadius", "rstar:kstar:mT:mult:part1:part2:mother1:mother2");

    for(auto iEv{0}; iEv<nEvents; ++iEv) {
        pythia.next();

        // evaluate multiplicity
        int nCh{0};
        for(auto iPart{0}; iPart<pythia.event.size(); ++iPart) {
            int pdg = std::abs(pythia.event[iPart].id());
            int status = std::abs(pythia.event[iPart].status());
            int eta = pythia.event[iPart].eta();
            if(IsStable(pdg) && std::abs(eta)<1 && ((status > 80 && status < 90) || (status > 111 && status < 119))) { // to be improved
                nCh++;
            }
        }

        for(auto iPart{2}; iPart<pythia.event.size(); ++iPart) {

            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);
            if(absPdg < 211 || absPdg > 2212)
                continue;

            // check source only for primary particles
            int idxMother = pythia.event[iPart].mother1();
            int pdgMother = pythia.event[idxMother].id();
            if(!IsPrimary(std::abs(pdgMother)))
                continue;
 
            int eta = pythia.event[iPart].eta();
            if (std::abs(eta)>0.5)
                continue;

            ROOT::Math::PxPyPzMVector mom(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
            ROOT::Math::XYZTVector prodvtx(pythia.event[iPart].xProd()*MM2FM, pythia.event[iPart].yProd()*MM2FM, pythia.event[iPart].zProd()*MM2FM, pythia.event[iPart].tProd()*MM2FM);
            if(absPdg == 211 || absPdg == 321 || absPdg == 2212) {
                pdgLight.push_back(pdg);
                momLight.push_back(mom);
                prodvtxLight.push_back(prodvtx);
                idxMotherLight.push_back(idxMother);
                pdgMotherLight.push_back(pdgMother);
            }
            else if(absPdg == 411 || absPdg == 421 || absPdg == 413) {
                pdgCharm.push_back(pdg);
                momCharm.push_back(mom);
                prodvtxCharm.push_back(prodvtx);
                idxMotherCharm.push_back(idxMother);
                pdgMotherCharm.push_back(pdgMother);
            }
            else {
                continue;
            }
        }

        // compute pairs of charm-light particles
        for(auto iCharm{0u}; iCharm<pdgCharm.size(); ++iCharm) {
            for(auto iLight{0u}; iLight<pdgLight.size(); ++iLight) {
                if(idxMotherCharm[iCharm] == idxMotherLight[iLight])
                    continue;
                float kStar{0}, rStar{0}, mT{0};
                ComputeRandKstar(momLight[iLight], momCharm[iCharm], prodvtxLight[iLight], prodvtxCharm[iCharm], kStar, rStar, mT);
                if(mT < 1.8 || kStar > 1)
                    continue;
                float arr4tuple[8] = {rStar, kStar, mT, (float)nCh, (float)pdgCharm[iCharm], (float)pdgLight[iLight], (float)pdgMotherCharm[iCharm], (float)pdgMotherLight[iLight]};
                tupleRadius->Fill(arr4tuple);
            }
        }

        for(auto iLightFirst{0u}; iLightFirst<pdgLight.size(); ++iLightFirst) {
            for(auto iLightSecond{iLightFirst+1}; iLightSecond<pdgLight.size(); ++iLightSecond) {
                if(idxMotherLight[iLightFirst] == idxMotherLight[iLightSecond])
                    continue;
                float kStar{0}, rStar{0}, mT{0};
                ComputeRandKstar(momLight[iLightFirst], momLight[iLightSecond], prodvtxLight[iLightFirst], prodvtxLight[iLightSecond], kStar, rStar, mT);
                if(mT < 1.8 || kStar > 1.)
                    continue;
                float arr4tuple[8] = {rStar, kStar, mT, (float)nCh, (float)pdgLight[iLightSecond], (float)pdgLight[iLightFirst], (float)pdgMotherLight[iLightSecond], (float)pdgMotherLight[iLightFirst]};
                tupleRadius->Fill(arr4tuple);
            }
        }

        pdgCharm.clear();
        pdgLight.clear();
        momCharm.clear();
        momLight.clear();
        prodvtxCharm.clear();
        prodvtxLight.clear();
        pdgMotherCharm.clear();
        pdgMotherLight.clear();
        idxMotherCharm.clear();
        idxMotherLight.clear();
    }

    TFile outFile(outFileName.data(), "recreate");
    tupleRadius->Write();
    outFile.Close();
}

//__________________________________________________________________________________________________
void ComputeRandKstar(ROOT::Math::PxPyPzMVector mom1, ROOT::Math::PxPyPzMVector mom2, ROOT::Math::XYZTVector vtx1, ROOT::Math::XYZTVector vtx2, float &kStar, float &rStar, float &mT) {
    ROOT::Math::PxPyPzMVector trackSum = mom1 + mom2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};

    ROOT::Math::PxPyPzMVector mom1CM = boostv12(mom1);
    ROOT::Math::PxPyPzMVector mom2CM = boostv12(mom2);

    ROOT::Math::XYZTVector vtx1CM = boostv12(vtx1);
    ROOT::Math::XYZTVector vtx2CM = boostv12(vtx2);

    ROOT::Math::PxPyPzMVector partRelK = mom1CM - mom2CM;
    kStar = 0.5 * partRelK.P();
    mT = std::sqrt(partRelK.M()*partRelK.M() + partRelK.Pt()*partRelK.Pt());

    ROOT::Math::XYZTVector partRelR = vtx1CM - vtx2CM;
    rStar = 0.5 * std::sqrt(partRelR.X()*partRelR.X() + partRelR.Y()*partRelR.Y() + partRelR.Z()*partRelR.Z());
}

//__________________________________________________________________________________________________
bool IsStable(int absPdg) {
    if(absPdg == 11 || absPdg == 13 || absPdg == 211 || absPdg == 321  || absPdg == 2212 || absPdg > 1000000000)
        return true;

    return false;
}

//__________________________________________________________________________________________________
bool IsPrimary(int absPdg) {
    if(absPdg >= 1 && absPdg <= 6)
        return true;
    if(absPdg == 0 || absPdg == 21)
        return true;
    if(absPdg == 1103 || absPdg == 2101 || absPdg == 2103 || absPdg == 2203)
        return true;
    if(absPdg == 3101 || absPdg == 3103 || absPdg == 3201 || absPdg == 3203 || absPdg == 3303)
        return true;
    if(absPdg == 4101 || absPdg == 4103 || absPdg == 4201 || absPdg == 4203 || absPdg == 4301 || absPdg == 4303 || absPdg == 4403)
        return true;
    if(absPdg == 5101 || absPdg == 5103 || absPdg == 5201 || absPdg == 5203 || absPdg == 5301 || absPdg == 5303 || absPdg == 5401 || absPdg == 5403 || absPdg == 5503)
        return true;

    return false;
}

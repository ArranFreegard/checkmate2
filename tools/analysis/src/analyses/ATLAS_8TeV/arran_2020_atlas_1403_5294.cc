#include "arran_2020_atlas_1403_5294.h"
// AUTHOR: Arran Freegard
//  EMAIL: acf1g14@soton.ac.uk
void Arran_2020_atlas_1403_5294::initialize() {
  setAnalysisName("arran_2020_atlas_1403_5294");          
  setInformation(""
    "# direct production of charginos, neutralinos and sleptons\n"
    "# two leptons and missing transverse momentum\n"
    "# 8 TeV, 20.3 fb^-1\n"
  "");
  setLuminosity(20.3*units::INVFB);      
  bookSignalRegions("mT2_90_SF;mT2_90_DF;mT2_120_SF;mT2_120_DF;mT2_150_SF;mT2_150_DF;WWa_SF;WWa_DF;WWb_SF;WWb_DF;WWc_SF;WWc_DF;Zjets");

}

void Arran_2020_atlas_1403_5294::analyze() {

  missingET->addMuons(muonsCombined);  // Adds muons to missing ET. This should almost always be done which is why this line is not commented out.

//Signal electrons fit the 'tight' criteria (defined in paper)
//Electrons with pT>10 GeV, |eta|<2.47 (defined in paper)
  electronsTight = filterPhaseSpace(electronsTight, 10., -2.47, 2.47);

//Muons with pT>10 GeV, |eta|<2.4 (defined in paper)
  muonsCombined = filterPhaseSpace(muonsCombined, 10., -2.4, 2.4);  

//Jets with pT>20 GeV, |eta|<4.5 (defined in paper)
  jets = filterPhaseSpace(jets, 20., -4.5, 4.5); 
 
//Standard isolation condition filters but not used, do my own isolation instead later
  //electronsLoose = filterIsolation(electronsLoose);
  //electronsTight = filterIsolation(electronsTight);
  //muonsCombined = filterIsolation(muonsCombined);


  std::vector<Jet*> bjets;
  std::vector<Jet*> lightjets;

//  std::vector<Jet*> overlapRemoval_muon_jet_tracks(std::vector<Jet*> cand_jets, std::vector<Muon*> cand_muons, double deltaR, int nTracks);


  // 1) If two electrons overlap within deltaR < 0.05, only the harder electron is stored (defined in paper, pg7)
electronsTight = overlapRemoval(electronsTight, 0.05, "y");

  // 2) Removes all jets for which there exists any electron with deltaR < 0.2. (defined in paper, pg7)
jets = overlapRemoval(jets, electronsTight, 0.2, "y");

//3) Remove taus if electron/muons within deltaR<0.2 (defined in paper, pg7)
//taujets = overlapRemoval(taujet, electronsLoose, 0.2, "y");
//taujets = overlapRemoval(taujet, electronsTight, 0.2, "y");
//taujets = overlapRemoval(taujet, muonsCombined, 0.2, "y");

  //4) Removes all electrons for which there exists any jet with deltaR < 0.4 (defined in paper, pg7)
electronsTight = overlapRemoval(electronsTight, jets, 0.4, "y");

  //5) Removes all muons for which there exists any jet with deltaR < 0.4 (defined in paper, pg7)
muonsCombined = overlapRemoval(muonsCombined, jets, 0.4, "y");   

  //6) Removes all electrons for which there exists any muon with deltaR < 0.01 (defined in paper, pg7)
electronsTight = overlapRemoval(electronsTight, muonsCombined, 0.01, "y");   

  //7) Removes all muons for which there exists any electron with deltaR < 0.01 (defined in paper, pg7)
muonsCombined = overlapRemoval(muonsCombined, electronsLoose, 0.01, "y");  
muonsCombined = overlapRemoval(muonsCombined, electronsTight, 0.01, "y"); 

//Remove both muons for deltaR < 0.05 (defined in paper, pg7)
//muonsCombined = overlapRemoval(muonsCombined, muonsCombined, 0.05);
muonsCombined = overlapRemoval(muonsCombined, 0.05);

//mll>20 is imposed later, no need for mll>12GeV

//Remove jets if tau within deltaR<0.2 (defined in paper, pg7)
//lightjets = overlapRemoval(lightjets, taujet, 0.2, "y");
//bjets = overlapRemoval(bjets, taujet, 0.2, "y");


//bjets with |eta|<2.4, or light jet with |eta|<2.4
    for (int i = 0; i < jets.size(); i++) 
    if ( fabs(jets[i]->Eta) < 2.4 && checkBTag(jets[i]) ) bjets.push_back(jets[i]);
    else lightjets.push_back(jets[i]);  


//============================================================================================================================================================
//============================================================================================================================================================
  
  std::vector<Electron*> electrons_base = electronsLoose;
  std::vector<Muon*> muons_base = muonsCombined;
    
//Isolation requirements are applied to signal electrons:
//The scalar sum of the p_T of tracks within a variable-sized cone around the lepton, excluding its own track, must be less than 16% of the lepton p_T. 
  //The track isolation cone radius for electrons is given by the smaller of ∆R = 10 GeV/p_T and ∆R = 0.2, that is,
  //a cone of size 0.2 at low p T but narrower for high-p T leptons. 
  //The energy of calorimeter energy clusters in a cone of ∆R = 0.2 around the electron (excluding the deposition
  //from the electron itself) must be less than 18% of the electron p T

//Electron isolation criteria (defined in paper, pg7)
  std::vector<Electron*> electrons_signal = Isolate_leptons_with_inverse_track_isolation_cone(electronsTight, tracks, towers, 0.3, 10, 0.3, 0.16, 0.18, true);


  //The track isolation cone radius for muons is given by the smaller of ∆R = 10 GeV/p_T and ∆R = 0.3, that is,
  //a cone of size 0.3 at low p T but narrower for high-p T leptons. 

//Muon isolation criteria  (defined in paper, pg7)
  std::vector<Muon*> muons_signal = Isolate_leptons_with_inverse_track_isolation_cone(muonsCombined, tracks, towers, 0.3, 10, 0.3, 0.16, 0.18, false);


  countCutflowEvent("00_all"); 

  double met = missingET->P4().Et();

//if missing energy is ~0, remove
  if ( met < 0.01) return;  
//  countCutflowEvent("01_MET>0"); 

//if the total number of electrons and muons is not 2, remove  
  if ( (muons_base.size() + electrons_base.size()) != 2 ) return;
  countCutflowEvent("01_2_base_leptons"); 

//Triggers, ee:97%, mumu:89%, emu:75% (defined in paper, pg 8)
  double triggerRatio = (double) rand() / (double) (RAND_MAX + 1.);

//If there is two electrons, but random number is greater than 97%, remove
  if ( electrons_signal.size() == 2 && triggerRatio > 0.97 ) return;

//If there is two muons, but random number is greater than 89%, remove
  if ( muons_signal.size() == 2 && triggerRatio > 0.89 ) return;

//If there is one electron and one muon, and random number is greater than 75%, remove
  if ( electrons_signal.size() == 1 && muons_signal.size() == 1 && triggerRatio > 0.75 ) return;

//if the total number of signal electrons and muons is not 2, remove    

//If leading lepton pT is less than 35 GeV or second lepton pT is less than 20 GeV, remove
//if ( muons_signal[0]->P4().Perp() < 35. or muons_signal[1]->P4().Perp() < 20. and electrons_signal[0]->P4().Perp() < 35. or electrons_signal[1]->P4().Perp() < 20. and muons_signal[0]->P4().Perp() < 35. or electrons_signal[1]->P4().Perp() < 20. and electrons_signal[0]->P4().Perp() < 35. or muons_signal[1]->P4().Perp() < 20.) return;

  bool El = false; bool Mu = false; bool ElMu = false;
  std::string flavour;
  std::vector<TLorentzVector> leptons;
  std::string tag;
//  bool sc = false;


//If the total electron + muon number is not 2, remove
  if (electrons_signal.size() + muons_signal.size() != 2) return;
//If the total electron number is 2, label E and add two leading electrons
  if (electrons_signal.size() == 2) {
    El = true;
    flavour = "El";
    leptons.push_back( electrons_signal[0]->P4() );
    leptons.push_back( electrons_signal[1]->P4() );
  }
//If the total muon number is 2, label M and add two leading muons
  else if (muons_signal.size() == 2) {
    Mu = true;
    flavour = "Mu";
    leptons.push_back( muons_signal[0]->P4() );
    leptons.push_back( muons_signal[1]->P4() );
  }
//If neither of these two cases, label EM and add leading muon and electron
  else {
    ElMu = true;
    flavour = "ElMu";
//If Electron is high pT than muon, add electron as leading
    if ( electrons_signal[0]->PT > muons_signal[0]->PT ) {
      leptons.push_back( electrons_signal[0]->P4() );
      leptons.push_back( muons_signal[0]->P4() );
    }
//Otherwise add muon as leading
    else {
      leptons.push_back( muons_signal[0]->P4() );
      leptons.push_back( electrons_signal[0]->P4() );
    }
  }

  countCutflowEvent("02_2_signal_leptons"); 


//  countCutflowEvent("05_opposite_charge");
  
//If leading lepton pT is less than 35 GeV or second lepton pT is less than 20 GeV, remove
  if ( leptons[0].Perp() < 35. or leptons[1].Perp() < 20.) return;  
  countCutflowEvent("03_2OS_lep1_pt>35,20");

//Remove if mll<20 GeV
  double mll = (leptons[0] + leptons[1]).M();
  if ( mll < 20.) return;
  countCutflowEvent("04_mll>20"); 

//  double mtau = mtautau(leptons[0], leptons[1], missingET->P4());
//  if ( mtau >= 0. ) return;
//  countCutflowEvent("05_tauveto");


//Remove if the jet size is not zero
  if (jets.size() != 0 ) return;
  countCutflowEvent("06_jetveto");

//If |mll-mZ| is less than 10 GeV, remove these
  if (fabs( mll-91.118 ) < 10. ) return;  
  countCutflowEvent("07_Zveto");

//WWa_Conversion (mll<120, EtMissRel>80, ptll > 80) 
//If mll>120, remove these 
// if ( mll > 120.) { return;
      //countCutflowEvent("08_SRWWa_mll_120");
//If pTll<80, remove theses
//     if ((leptons[0] + leptons[1]).Perp() < 80.) {return;
          // countCutflowEvent("09_SRWWa_ptll>80");
//If MET<80, remove theses
//      if (missingET->P4().Et() < 80.) {return;
//             countCutflowEvent("10_SRWWa_ETMissRel>80"); 
//             countSignalEvent("WWa_SF");
//      }
 //   }
 //     }
if (mll < 120. and (leptons[0] + leptons[1]).Perp() > 80. and missingET->P4().Et() > 80.) {
   if (El or Mu ) countSignalEvent("WWa_SF");
   else countSignalEvent("WWa_DF");
   }



//WWb
//If mT2 is less than 90, remove these
//  if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) < 90. ) {
//     countCutflowEvent("11_SRWWb_mt2_90");  

//If mll>170, remove these events   
//     if ( mll > 170. ) {
//     	countCutflowEvent("12_SRWWb_mll_170");
//        countSignalEvent("WWb_SF");
//                       }
//                                                }
if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90.  and  mll < 170. ) {
    if (El or Mu ) countSignalEvent("WWb_SF");
    else countSignalEvent("WWb_DF");
    }

//WWc
//If mT2 is less than 100, remove these
//  if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) < 100. )  {
//       countCutflowEvent("13_SRWWc_mt2_100");
//       countSignalEvent("WWc_SF");

                                     
if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 100. ) {
   if (El or Mu ) countSignalEvent("WWc_SF");
   else countSignalEvent("WWc_DF");
   }

//mT2_90
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent(flavour+"_mT2_90");
      if (El or Mu ) countSignalEvent("mT2_90_SF");
      else countSignalEvent("mT2_90_DF");
    }

//mt2_120
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 120. ) {
      countCutflowEvent(flavour+"_mT2_120");
      if (El or Mu ) countSignalEvent("mT2_120_SF");
      else countSignalEvent("mT2_120_DF");
    }
 
//mt2_150
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 150. ) {
      countCutflowEvent(flavour+"_mT2_150");
      if (El or Mu ) countSignalEvent("mT2_150_SF");
      else countSignalEvent("mT2_150_DF");
   }

//zjet
//    if ( btag || forwardjet ) return;
//    countCutflowEvent(flavour+"_b_veto");
// "Z window" line in Zjets SRs
//    if ( fabs(mll - 91.2) > 10. ) return;
//    countCutflowEvent(flavour+"_Z_window");    
//    if ( pll.Perp() < 80. ) return;
//    countCutflowEvent(flavour+"_pTll");  
//    if (missingETrel < 80. ) return;
//    countCutflowEvent(flavour+"_Etmiss");
// "0.3 < \Delta R_{ll} < 1.5" in Zjets SRs
//    if ( leptons[0].DeltaR(leptons[1]) < 0.3 || leptons[0].DeltaR(leptons[1]) > 1.5 ) return;
//    countCutflowEvent(flavour+"_DRll");
// "50 < m_{jj} < 100 GeV" in Zjets SRs    
//    double mjj = (centraljets[0] + centraljets[1]).M();
//    if ( mjj < 50. || mjj > 100. ) return;
//    countCutflowEvent(flavour+"_mjj");
    // "Jet p_T > 45 GeV" in Zjets SRs corresponds to "The two highest-pT central light jets must have pT > 45 GeV"    
//    if ( centraljets[1].Perp() < 45. ) return;
//    countCutflowEvent(flavour+"_jet_pT");  
//    countSignalEvent("Zjets");
//  }
//  else return;
  }

void Arran_2020_atlas_1403_5294::finalize() {

}       

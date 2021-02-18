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

  std::vector<Jet*> bjets;
  std::vector<Jet*> lightjets;
  std::vector<Jet*> forwardjets;

///OVERLAP REMOVALS///
  electronsTight = overlapRemoval(electronsTight, 0.05, "y");
  jets = overlapRemoval(jets, electronsTight, 0.2, "y");
  electronsTight = overlapRemoval(electronsTight, jets, 0.4, "y");
  muonsCombined = overlapRemoval(muonsCombined, jets, 0.4, "y");   

//muonsCombined = overlapRemoval(muonsCombined, electronsTight, 0.01, "y"); 
//muonsCombined = overlapRemoval(muonsCombined, 0.05, "y");


///JETS///
//bjets with |eta|<2.4, or light jet with |eta|<2.4
  for (int i = 0; i < jets.size(); i++) {
    if ( (fabs(jets[i]->Eta) < 2.4) && (jets[i]->PT > 20.) ){
      if (checkBTag(jets[i]) ) bjets.push_back(jets[i]);
      else {
        lightjets.push_back(jets[i]);
        }
      }
    else if (jets[i]->PT > 30. )  forwardjets.push_back(jets[i]);
    }

///ISOLATION REQUIREMENTS///
  std::vector<Electron*> electrons_base = electronsTight;
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

////////////////////////////////////////////////////////////
  countCutflowEvent("01_Initial"); 

  bool ee = false; bool mumu = false; bool emu = false;
  std::string flavour;
  std::vector<TLorentzVector> leptons;
  if (electrons_signal.size() + muons_signal.size() != 2) return;
  if (electrons_signal.size() == 2) {
    ee = true;
    flavour = "ee";
    leptons.push_back( electrons_signal[0]->P4() );
    leptons.push_back( electrons_signal[1]->P4() );
  }
  else if (muons_signal.size() == 2) {
    mumu = true;
    flavour = "mumu";
    leptons.push_back( muons_signal[0]->P4() );
    leptons.push_back( muons_signal[1]->P4() );
  }
  else {
    emu = true;
    flavour = "emu";
    if ( electrons_signal[0]->PT > muons_signal[0]->PT ) {
      leptons.push_back( electrons_signal[0]->P4() );
      leptons.push_back( muons_signal[0]->P4() );
    }
    else {
      leptons.push_back( muons_signal[0]->P4() );
      leptons.push_back( electrons_signal[0]->P4() );
    }
  }
  
  TLorentzVector pll = leptons[0] + leptons[1];
  double mll = pll.M();

  if ( leptons[0].Perp() < 35. || leptons[1].Perp() < 20. ) return;

  countCutflowEvent("02_2_OS_leptons"); 

  if ( mll < 20.) return;

  countCutflowEvent("03_mll>20_GeV");

  countCutflowEvent("04_tau_veto");

//ee-channel//
//////////////
  if ( ee ) {
    countCutflowEvent("05_ee_leptons");
    if (jets.size() != 0 ) return;
    countCutflowEvent("06_ee_jet_veto");
    if (fabs( mll-91.118 ) < 10. ) return;  
    countCutflowEvent("07_ee_Z_veto");
//Wwa
    if ((leptons[0] + leptons[1]).Perp() > 80.) {
      countCutflowEvent("08_ee_WWa_pTll>80_GeV");
      if (missingET->P4().Et() > 80.)
        countCutflowEvent("09_ee_WWa_METrel>80_GeV");
        if (mll < 120.)
          countCutflowEvent("10_ee_WWa_mll<120_GeV");
          countSignalEvent("WWa_SF");
    }

//Wwb
    if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent("08_ee_WWb_mT2>90_GeV");
      if (mll < 170.)
      countCutflowEvent("09_ee_WWb_mll<170_GeV");
      countSignalEvent("WWb_SF");
    }
//Wwc
    if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 100. ) {
      countCutflowEvent("08_ee_WWc_mT2>100_GeV");
      countSignalEvent("WWc_SF");
    }

//mT2_90
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent("08_ee_mT2>90_GeV");
      countSignalEvent("mT2_90_SF");
    }

//mT2_120
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 120. ) {
      countCutflowEvent("08_ee_mT2>120_GeV");
      countSignalEvent("mT2_120_SF");
    }
//mT2_150
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 150. ) {
      countCutflowEvent("08_ee_mT2>150_GeV");
      countSignalEvent("mT2_150_SF");
    }
//Zjets

    if ( lightjets.size() > 1 ) {
      countCutflowEvent("07_ee_Zjets_>=2_central_light_jets");
    if ( bjets.size() != 0 || forwardjets.size() != 0 ) return;
      countCutflowEvent("08_ee_Zjets_b_&_forward_jet_veto");
    if ( fabs(mll - 91.2) > 10. ) return;
    countCutflowEvent("09_ee_Zjets_Z_window");
    if ( pll.Perp() < 80. ) return;
    countCutflowEvent("10_ee_Zjets_pTll>80_GeV");
    if ( missingET->P4().Et() < 80. ) return;
    countCutflowEvent("11_ee_Zjets_MET_rel>80_GeV");
    if ( leptons[0].DeltaR(leptons[1]) < 0.3 || leptons[0].DeltaR(leptons[1]) > 1.5 ) return;
    countCutflowEvent("12_ee_Zjets_0.3<dRll<1.5");
    double mjj = (lightjets[0]->P4() + lightjets[1]->P4()).M();
    if ( mjj < 50. || mjj > 100. ) return;
    countCutflowEvent("13_ee_Zjets_50<mjj<100");
    if ( lightjets[1]->P4().Perp() < 45. ) return;
    countCutflowEvent("14_ee_Zjets_pTjj>45");
    countSignalEvent("Zjets");
    }
  }

//mumu-channel//
////////////////
  else if ( mumu ) {
    countCutflowEvent("05_mumu_leptons");
    if (jets.size() != 0 ) return;
    countCutflowEvent("06_mumu_jet_veto");
    if (fabs( mll-91.118 ) < 10. ) return;  
    countCutflowEvent("07_mumu_Z_veto");
//Wwa
    if ((leptons[0] + leptons[1]).Perp() > 80.) {
      countCutflowEvent("08_mumu_WWa_pTll>80_GeV");
      if (missingET->P4().Et() > 80.)
        countCutflowEvent("09_mumu_WWa_METrel>80_GeV");
        if (mll < 120.)
          countCutflowEvent("10_mumu_WWa_mll<120_GeV");
          countSignalEvent("WWa_SF");
    }

//Wwb
    if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent("08_mumu_WWb_mT2>90_GeV");
      if (mll < 170.)
      countCutflowEvent("09_mumu_WWb_mll<170_GeV");
      countSignalEvent("WWb_SF");
    }
//Wwc
    if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 100. ) {
      countCutflowEvent("08_mumu_WWc_mT2>100_GeV");
      countSignalEvent("WWc_SF");
    }

//mT2_90
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent("08_mumu_mT2>90_GeV");
      countSignalEvent("mT2_90_SF");
    }

//mT2_120
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 120. ) {
      countCutflowEvent("08_mumu_mT2>120_GeV");
      countSignalEvent("mT2_120_SF");
    }
//mT2_150
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 150. ) {
      countCutflowEvent("08_mumu_mT2>150_GeV");
      countSignalEvent("mT2_150_SF");
    }
//Zjets
    if ( lightjets.size() > 1 ) {
      countCutflowEvent("07_mumu_Zjets_>=2_central_light_jets");
    if ( bjets.size() != 0 || forwardjets.size() != 0 ) return;
      countCutflowEvent("08_mumu_Zjets_b_&_forward_jet_veto");
    if ( fabs(mll - 91.2) > 10. ) return;
    countCutflowEvent("09_mumu_Zjets_Z_window");
    if ( pll.Perp() < 80. ) return;
    countCutflowEvent("10_mumu_Zjets_pTll>80_GeV");
    if ( missingET->P4().Et() < 80. ) return;
    countCutflowEvent("11_mumu_Zjets_MET_rel>80_GeV");
    if ( leptons[0].DeltaR(leptons[1]) < 0.3 || leptons[0].DeltaR(leptons[1]) > 1.5 ) return;
    countCutflowEvent("12_mumu_Zjets_0.3<dRll<1.5");
    double mjj = (lightjets[0]->P4() + lightjets[1]->P4()).M();
    if ( mjj < 50. || mjj > 100. ) return;
    countCutflowEvent("13_mumu_Zjets_50<mjj<100");
    if ( lightjets[1]->P4().Perp() < 45. ) return;
    countCutflowEvent("14_mumu_Zjets_pTjj>45");
    countSignalEvent("Zjets");
    }
  }

//emu-channel//
///////////////
  else {
    countCutflowEvent("05_emu_leptons");
    if (jets.size() != 0 ) return;
    countCutflowEvent("06_emu_jet_veto");
//Wwa
    if ((leptons[0] + leptons[1]).Perp() > 80.) {
      countCutflowEvent("08_emu_WWa_pTll>80_GeV");
      if (missingET->P4().Et() > 80.)
        countCutflowEvent("09_emu_WWa_METrel>80_GeV");
        if (mll < 120.)
          countCutflowEvent("10_emu_WWa_mll<120_GeV");
          countSignalEvent("WWa_SF");
    }

//Wwb
    if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent("08_emu_WWb_mT2>90_GeV");
      if (mll < 170.)
      countCutflowEvent("09_emu_WWb_mll<170_GeV");
      countSignalEvent("WWb_SF");
    }
//Wwc
    if ( mT2( leptons[0], leptons[1], 0., missingET->P4()) > 100. ) {
      countCutflowEvent("08_emu_WWc_mT2>100_GeV");
      countSignalEvent("WWc_SF");
    }

//mT2_90
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 90. ) {
      countCutflowEvent("08_emu_mT2>90_GeV");
      countSignalEvent("mT2_90_SF");
    }

//mT2_120
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 120. ) {
      countCutflowEvent("08_emu_mT2>120_GeV");
      countSignalEvent("mT2_120_SF");
    }
//mT2_150
    if (mT2( leptons[0], leptons[1], 0., missingET->P4()) > 150. ) {
      countCutflowEvent("08_emu_mT2>150_GeV");
      countSignalEvent("mT2_150_SF");
    }
  }
}

void Arran_2020_atlas_1403_5294::finalize() {

} 

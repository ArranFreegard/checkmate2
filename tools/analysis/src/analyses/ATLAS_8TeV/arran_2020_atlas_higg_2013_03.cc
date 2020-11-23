#include "arran_2020_atlas_higg_2013_03.h"
// AUTHOR: Arran Freegard
//  EMAIL: acf1g14@soton.ac.uk
void Arran_2020_atlas_higg_2013_03::initialize() {
  setAnalysisName("arran_2020_atlas_higg_2013_03");          
  setInformation(""
    "# no\n"
  "");
  setLuminosity(20.3*units::INVFB);      
  bookSignalRegions("sig");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}
 // All std::vector members and etmiss have the common properties PT, Eta, Phi and P4() with the latter giving access to the full ROOT TLorentzVector.
void Arran_2020_atlas_higg_2013_03::analyze() {
//Signal electrons fit the 'tight' criteria (defined in paper)
//Electrons with pT>7 GeV, |eta|<2.47 (defined in paper)
  electronsLoose = filterPhaseSpace(electronsLoose, 7., -2.47, 2.47);
  electronsTight = filterPhaseSpace(electronsTight, 7., -2.47, 2.47);

//Muons with pT>7 GeV, |eta|<2.5 (defined in paper)
  muonsCombined = filterPhaseSpace(muonsCombined, 7., -2.4, 2.4);  

//Jets with pT>25 GeV, |eta|<2.5 (defined in paper)
  jets = filterPhaseSpace(jets, 25., -2.5, 2.5); 
 

  tracks = filterPhaseSpace(tracks, 0.5, -2.5, 2.5);

//discriminate against jets from additional minimum bias interactions, selection criteria
//are applied to ensure that most of the jet momentum, for jets with |η| < 2.5, is 
//associated with tracks originating from the primary vertex, which is taken to be the 
//vertex with the highest summed p^2T of associated tracks.


//Standard isolation condition filters but not used, do my own isolation instead later
  //electronsLoose = filterIsolation(electronsLoose);
  //electronsTight = filterIsolation(electronsTight);
  //muonsCombined = filterIsolation(muonsCombined);


  //1) Removes all electrons for which there exists any muon with deltaR < 0.2 (defined in paper, pg2)
electronsLoose = overlapRemoval(electronsLoose, muonsCombined, 0.2, "y");  
electronsTight = overlapRemoval(electronsTight, muonsCombined, 0.2, "y"); 

  // 2) Removes all jets for which there exists any electron with deltaR < 0.2. (defined in paper, pg2)
jets = overlapRemoval(jets, electronsLoose, 0.2, "y");
jets = overlapRemoval(jets, electronsTight, 0.2, "y");

  //3) Removes all electrons for which there exists any jet with deltaR < 0.4 (defined in paper, pg2)
electronsLoose = overlapRemoval(electronsLoose, jets, 0.4, "y");
electronsTight = overlapRemoval(electronsTight, jets, 0.4, "y");

  //4) Removes all muons for which there exists any jet with deltaR < 0.4 (defined in paper, pg2)
muonsCombined = overlapRemoval(muonsCombined, jets, 0.4, "y");   

//============================================================================================================================================================
//============================================================================================================================================================

  missingET->addMuons(muonsCombined);  // Adds muons to missing ET. This should almost always be done which is why this line is not commented out.
  
  std::vector<Electron*> electrons_base = electronsLoose;
  std::vector<Muon*> muons_base = muonsCombined;
    
//Isolation requirements are applied to signal electrons:
//The scalar sum of the p_T of tracks within a variable-sized cone around the lepton, excluding its own track, must be less than 10% of the lepton p_T. 
  //The track isolation cone radius for electrons is given by the smaller of ∆R = 10 GeV/p_T and ∆R = 0.2, that is,
  //a cone of size 0.2 at low p T but narrower for high-p T leptons. 
  //The energy of calorimeter energy clusters in a cone of ∆R = 0.2 around the electron (excluding the deposition
  //from the electron itself) must be less than 18% of the electron p T

//Electron isolation criteria (defined in paper, pg7)
  std::vector<Electron*> electrons_signal = Isolate_leptons_with_inverse_track_isolation_cone(electronsTight, tracks, towers, 0.2, 10, 0.2, 0.10, 0.10, false);


  //The track isolation cone radius for muons is given by the smaller of ∆R = 10 GeV/p_T and ∆R = 0.2, that is,
  //a cone of size 0.2 at low p T but narrower for high-p T leptons. 

//Muon isolation criteria  (defined in paper, pg7)
  std::vector<Muon*> muons_signal = Isolate_leptons_with_inverse_track_isolation_cone(muonsCombined, tracks, towers, 0.2, 10, 0.2, 0.10, 0.10, false);

//  countCutflowEvent("00_all"); 

  double met = missingET->P4().Et();


//if the total number of electrons and muons is not 2, remove  
  if ( (muons_base.size() + electrons_base.size()) != 2 ) return;
//  countCutflowEvent("01_2_base_leptons"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Unknown trigger efficiencies!
//Triggers, ee:97%, mumu:89%, emu:75% (defined in paper, pg 8)
//  double triggerRatio = (double) rand() / (double) (RAND_MAX + 1.);
//If there is two electrons, but random number is greater than 97%, remove
//  if ( electrons_signal.size() == 2 && triggerRatio > 0.97 ) return;
//If there is two muons, but random number is greater than 89%, remove
//  if ( muons_signal.size() == 2 && triggerRatio > 0.89 ) return;
//If there is one electron and one muon, and random number is greater than 75%, remove
//  if ( electrons_signal.size() == 1 && muons_signal.size() == 1 && triggerRatio > 0.75 ) return;
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//if the total number of signal electrons and muons is not 2, remove    
  if ( (muons_signal.size() + electrons_signal.size()) != 2 ) return;  


//  countCutflowEvent("01_2_signal_leptons"); 
  

  std::vector<TLorentzVector> leptons;
  std::string tag;
  bool sc = false;


//If, for two muons, the multiplied charge is -1, they are opposite charge
  if ( muons_signal.size() == 2 ) {
    leptons.push_back(muons_signal[0]->P4());
    leptons.push_back(muons_signal[1]->P4());
    tag = "M";
    if (muons_signal[0]->Charge * muons_signal[1]->Charge < 0) sc = true; 
  }

//If, for two electrons, the multiplied charge is -1, they are opposite charge
  else if ( electrons_signal.size() == 2 ) {
    leptons.push_back(electrons_signal[0]->P4());
    leptons.push_back(electrons_signal[1]->P4());
    tag = "E";
    if (electrons_signal[0]->Charge * electrons_signal[1]->Charge < 0) sc = true; 
  }  

  else if ( electrons_signal.size() == 1 && muons_signal.size() == 1 ) {
    leptons.push_back(electrons_signal[0]->P4());
    leptons.push_back(muons_signal[0]->P4());
    tag = "EMu";
    if (electrons_signal[0]->Charge * muons_signal[0]->Charge < 0) sc = true; 
  }  
//Otherwise, remove them
  else return;
 
//  countCutflowEvent("02_same_flav");

//If events do not fit the opposite charge tag, remove them 
  if ( !sc ) return;  
//  countCutflowEvent("05_opposite_charge");
  
//If leading lepton pT is less than 20 GeV, remove (pg2)
  if ( leptons[0].Perp() < 20. or leptons[1].Perp() < 20.) return;  
  countCutflowEvent("01_2OS_lep1_pt>20");
////////////////////////////////////////////////////////////////////////

  double mll = (leptons[0] + leptons[1]).M();

//"Z window");
//need 76 GeV < mll < 106 GeV
if (mll < 76. or mll > 106) return;  
   countCutflowEvent("02_Zwindow");

//"MET > 90 GeV");
if ( met < 90) return;
  countCutflowEvent("03_MET>90 GeV");

////////////////////////////////////////////////////////////////////
    // basic variables for the two leptons
    //TLorentzVector p_lepton1 = leptons[0].P();
    //TLorentzVector p_lepton2 = leptons[1].P();
    TLorentzVector p_dilepton = (leptons[0] + leptons[1]);
    //double mll = p_dilepton.M();
	//TLorentzVector pTll = (leptons[0] + leptons[1]).Perp();
    //TLorentzVector pTll = p_dilepton.Pt();

    // MET and pTmiss
    //TLorentzVector p_met = event.rec()->MET().momentum();
    //double met = p_met.Pt();
    //TLorentzVector ptmiss;
    //for(unsigned int ii=0; ii<myTracks.size(); ii++)
    //{
    //  ptmiss -= myTracks[ii]->momentum();

  TLorentzVector track_MET;
  for (int i = 0; i < tracks.size(); i++) 
     if ( tracks[i]->PT > 0.5 and fabs(tracks[i]->Eta) < 2.5 ) track_MET -= tracks[i]->P4();   
     //another calculation of MET 
  TLorentzVector pTmiss = TLorentzVector(0., 0., 0., 0.);
  for (int i = 0; i < jets.size(); i++)
    pTmiss -= jets[i]->P4();
  
  TLorentzVector pTmiss_soft = TLorentzVector(0., 0., 0., 0.);
  for (std::vector<Track*>::iterator it=tracks.begin(); it!=tracks.end(); it++) {
    bool used = false;
    for (int i = 0; i < jets.size() && !used; i++)
      for (int part = 0; part < jets[i]->Particles.GetEntries(); part++)
        if ((*it)->PT > 0.4 && jets[i]->Particles.At(part) == (*it)->Particle) {
	  used = true;
	  break;
	}
    if (!used)  pTmiss_soft -= (*it)->P4();
  }
  
  pTmiss += pTmiss_soft;
  
///////////////////////////////////////////////////////////////////

//"dilepton-MET separation");
    // azimuthal separation between...

	//double dphi_met_pTmiss = fabs(met.DeltaPhi(pTmiss)); // MET and pTmiss

//    double dphi_met_ptmiss = fabs(p_met.DeltaPhi(ptmiss)); // MET and pTmiss

    //double dphi_ll_met = fabs((leptons[0] + leptons[1]).DeltaPhi(met)); // di-lepton pair and MET

//    double dphi_ll_met = fabs(p_dilepton.DeltaPhi(p_met)); // di-lepton pair and MET
//    double dphi_l_l = fabs(p_lepton1.DeltaPhi(p_lepton2)); // the two leptons

    //double dphi_l_l = fabs(leptons[0].DeltaPhi(leptons[1])); // the two leptons

// similarity between pT of dilepton pair and MET
//    double frac_diff = fabs(met - pTll) / pTll;

//	double frac_diff = fabs(  ( track_MET - (leptons[0] + leptons[1]).Perp() ) / ( (leptons[0] + leptons[1]).Perp() )  );

//"pTmiss-MET separation");
    // cut on the azimuthal angle between MET and pTmiss to reject MET from
    // misreconstructed energy in the calorimeter (dphi < 0.2)
	//if (dphi_met_pTmiss > 0.2) return;    
	if (fabs(track_MET.DeltaPhi(pTmiss)) > 0.2) return;
 countCutflowEvent("04_pTmiss-MET_separation");

//"dilepton-MET separation"
//    the azimuthal separation between the dilepton system and the MET
    // is required to be larger than 2.6

    //if (dphi_ll_met < 2.6) return;
    if (fabs((leptons[0] + leptons[1]).DeltaPhi(missingET->P4())) < 2.6) return;
 //   if (fabs((leptons[0] + leptons[1]).DeltaPhi(met)) < 2.6) return;
 countCutflowEvent("05_dilepton-MET_separation");

//"lepton-lepton separation");
    // the azimuthal opening angle of the two leptons is required to be less
    // than 1.7

    //if (dphi_l_l > 1.7)  return;
    if (fabs(leptons[0].DeltaPhi(leptons[1])) > 1.7)  return;
 countCutflowEvent("06_lepton-lepton_separation");


//"pTll-MET similarity");
    // pT of the dilepton pair and MET should be similar
    if (fabs(  (met - p_dilepton.Perp()) / p_dilepton.Perp() ) > 0.2) return;
 countCutflowEvent("07_pTll-MET_similarity");


//"jet veto");
//Remove if the jet size is not zero
  if (jets.size() != 0 ) return;
	  countCutflowEvent("08_jetveto");

  		countSignalEvent("sig");
/////////////////////////////////////////////////////

    }

void Arran_2020_atlas_higg_2013_03::finalize() {
  // Whatever should be done after the run goes here
}       

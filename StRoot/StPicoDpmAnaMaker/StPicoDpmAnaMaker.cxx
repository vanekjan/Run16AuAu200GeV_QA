#include "StPicoDpmAnaMaker.h"
//#include "StPicoHFMaker/StHFCuts.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TMath.h"

ClassImp(StPicoDpmAnaMaker)

using namespace std;

// _________________________________________________________
StPicoDpmAnaMaker::StPicoDpmAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  
				       char const* inputHFListHFtree = "") :
  StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
  mDecayChannel(kChannel1), mRefmultCorrUtil(NULL),mOutFileBaseName(outputBaseFileName){
   

  // constructor
}

// _________________________________________________________
StPicoDpmAnaMaker::~StPicoDpmAnaMaker() {
  // destructor
}

// _________________________________________________________
int StPicoDpmAnaMaker::InitHF() {
  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written

  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());

 mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

//------------------GET AND SAVE LIST OF RUNS------------------------------------------------------------

  ifstream RunList("./ListOfRuns/Runnumber.list");
  

  if(RunList.is_open())
  {

    int Run;
    while( RunList >> Run)
    {    
      RunNumberVector.push_back(Run);
    }
  }
  else
  {
    cout<<"Failed to open file!"<<endl;
  }


//----------------------------------------------------------------------------------------------


//-------MY Run16 QA HISTOGRAMS-----------------
	

//___detector, centrality and statistics (No. of events, tracks...) histograms
	mOutList->Add(new TH2F("h_mh1Cent", "EventsVsCentrality;cent;CountsvsRunIndex", 10, -1.5, 8.5, RunNumberVector.size()+1, -1, RunNumberVector.size()));
  mOutList->Add(new TH2F("h_mh1CentWg", "EventsVsCentrality;cent;CountsvsRunIndex", 10, -1.5, 8.5, RunNumberVector.size()+1, -1, RunNumberVector.size()));
  mOutList->Add(new TH2F("h_mh1gRefmultCor", "gRefmultCor;gRefmult;CountsvsRunIndex", 700, 0, 700, RunNumberVector.size()+1, -1, RunNumberVector.size()));
  mOutList->Add(new TH2F("h_mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;CountsvsRunIndex", 700, 0, 700, RunNumberVector.size()+1, -1, RunNumberVector.size()));
  mOutList->Add(new TH3F("h_mh2CentVz", "CentralityVsVz;cent;VzvsRunIndex", 10, -1.5, 8.5, 200, -10, 10, RunNumberVector.size()+1, -1, RunNumberVector.size()));
  mOutList->Add(new TH3F("h_mh2CentVzWg", "CentralityVsVzWg;cent;VzvsRunIndex", 10, -1.5, 8.5, 200, -10, 10, RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2D("h_QA_Vz", "Vz_vs_RunIndex", 200, -10, 10, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2D("h_QA_VzmVzVPD", "Vz-VzVPD_vs_RunIndex", 100, -4, 4, RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2D("h_QA_ZDC_rate", "ZDC_rateVsRunIndex", 200, 0, 100, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2D("h_QA_BBC_rate", "BBC_rateVsRunIndex", 200, 0, 100, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range

	mOutList->Add(new TH2D("h_QA_ZDC_rate_pileUp", "ZDC_rateVsRunIndex", 200, 0, 100, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2D("h_QA_BBC_rate_pileUp", "BBC_rateVsRunIndex", 200, 0, 100, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range

	mOutList->Add(new TH2D("h_QA_ZDC_rate_pileUp_TOF", "ZDC_rate_TOFVsRunIndex", 200, 0, 100, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2D("h_QA_BBC_rate_pileUp_TOF", "BBC_rate_TOFVsRunIndex", 200, 0, 100, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range

	mOutList->Add(new TH2I("h_QA_reweight_isNaN", "reweight_isNaNvsRunIndex", 2, 0, 2, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check for bad weights for refmutCorr
	mOutList->Add(new TH1I("h_QA_nEvents", "Number_of_eventsvsRunIndex", RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of events in run

	mOutList->Add(new TH2D("h_QA_nTracks", "Number_of_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size()));	//total nuber of tracks in event (no cuts)
	mOutList->Add(new TH2D("h_QA_nTracks_TPC", "Number_of_TPC_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of TPC tracks in run (heve to pass TPC cuts)
	mOutList->Add(new TH2D("h_QA_nTracks_HFT", "Number_of_HFT_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of HFT tracks (have to pass track->isHFTTrack())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_PXL1", "Number_of_HFT_PXL1_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of PXL1 tracks (have to pass track->hasPxl1Hit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_PXL2", "Number_of_HFT_PXL2_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of PXL2 tracks (have to pass track->hasPxl2Hit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_IST", "Number_of_HFT_IST_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size()));	//number of IST tracks (have to pass track->hasIstHit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_SSD", "Number_of_HFT_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size()));	//number of SSD tracks (have to pass track->hasSstHit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_IST_or_SSD", "Number_of_HFT_IST_or_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of IST or SST tracks (have to pass track->hasIstHit() || track->hasSstHit())
	mOutList->Add(new TH2D("h_QA_nTracks_TOF", "Number_of_TOF_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //nuber of TOF tracks (have to have TOF info)
	mOutList->Add(new TH2D("h_QA_nTracks_BEMC", "Number_of_BEMC_tracks_vs_RunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of BEMC tracks

	//_____________________HFT, eta cut_____________________________________________
	mOutList->Add(new TH2D("h_QA_nTracks_TPC_etaCut", "Number_of_TPC_etaCut_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut", "Number_of_HFT_etaCut_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of HFT tracks (have to pass track->isHFTTrack())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_PXL1", "Number_of_HFT_etaCut_PXL1_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of PXL1 tracks (have to pass track->hasPxl1Hit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_PXL2", "Number_of_HFT_etaCut_PXL2_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of PXL2 tracks (have to pass track->hasPxl2Hit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_IST", "Number_of_HFT_etaCut_IST_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size()));	//number of IST tracks (have to pass track->hasIstHit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_SSD", "Number_of_HFT_etaCut_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size()));	//number of SSD tracks (have to pass track->hasSstHit())
	mOutList->Add(new TH2D("h_QA_nTracks_HFT_etaCut_IST_or_SSD", "Number_of_HFT_etaCut_IST_or_SSD_tracksvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //number of IST or SST tracks (have to pass track->hasIstHit() || track->hasSstHit())

	//Number of tracks vs. pT - get from integral of pT spectrum

//____general QA histograms (all tracks within (TPC, track quality) cuts, NO TOF and HFT)__________________________________
	mOutList->Add(new TH2F("h_QA_pT", "Transverse_momentum_TPCvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2F("h_QA_eta", "Pseudorapidity_TPCvsRunIndex", 200, -1., 1., RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_phi", "Azimuthal_angle_TPCvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_vertex_x", "Vertex_x_positionvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_y", "Vertex_y_positionvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_z", "Vertex_z_positionvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_DCA_xy_TPC", "DCA_xy_TPCvsRunIndex", 600, -3, 3, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_TPC", "DCA_xy_zoom_TPCvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_z_TPC", "DCA_z_TPCvsRunIndex", 600, -3, 3, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_z_zoom_TPC", "DCA_z_zoom_TPCvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size()+1, -1, RunNumberVector.size()));


//___HFT QA hitograms______________________________________________________________________________________________________
	mOutList->Add(new TH2F("h_QA_pT_HFT", "Transverse_momentum_HFTvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2F("h_QA_eta_HFT", "Pseudorapidity_HFTvsRunIndex", 200, -1., 1., RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_phi_HFT", "Azimuthal_angle_HFTvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_vertex_x_HFT", "Vertex_x_position_HFTvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_y_HFT", "Vertex_y_position_HFTvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_z_HFT", "Vertex_z_position_HFTvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_DCA_xy_HFT", "DCA_xy_HFTvsRunIndex", 600, -3, 3, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_HFT", "DCA_xy_zoom_HFTvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_z_HFT", "DCA_z_HFTvsRunIndex", 600, -3, 3, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_z_zoom_HFT", "DCA_z_zoom_HFTvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size()+1, -1, RunNumberVector.size()));

//___HFT QA hitograms, eta cut______________________________________________________________________________________________________
	mOutList->Add(new TH2F("h_QA_pT_HFT_etaCut", "Transverse_momentum_HFT_etaCutvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2F("h_QA_eta_HFT_etaCut", "Pseudorapidity_HFT_etaCutvsRunIndex", 200, -1., 1., RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_phi_HFT_etaCut", "Azimuthal_angle_HFT_etaCutvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_vertex_x_HFT_etaCut", "Vertex_x_position_HFT_etaCutvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_y_HFT_etaCut", "Vertex_y_position_HFT_etaCutvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_z_HFT_etaCut", "Vertex_z_position_HFT_etaCutvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_DCA_xy_HFT_etaCut", "DCA_xy_HFT_etaCutvsRunIndex", 600, -3, 3, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_xy_zoom_HFT_etaCut", "DCA_xy_zoom_HFT_etaCutvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_z_HFT_etaCut", "DCA_z_HFT_etaCutvsRunIndex", 600, -3, 3, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_DCA_z_zoom_HFT_etaCut", "DCA_z_zoom_HFT_etaCutvsRunIndex", 100, -0.1, 0.1, RunNumberVector.size()+1, -1, RunNumberVector.size()));

//____TOF QA histograms______________________________________________________________________________________________________
	mOutList->Add(new TH2F("h_QA_pT_TOF", "Transverse_momentum_TOFvsRunIndex", 200, 0, 20, RunNumberVector.size()+1, -1, RunNumberVector.size())); //check binning and range
	mOutList->Add(new TH2F("h_QA_eta_TOF", "Pseudorapidity_TOFvsRunIndex", 200, -1., 1., RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_phi_TOF", "Azimuthal_angle_TOFvsRunIndex", 200, -TMath::Pi(), TMath::Pi(), RunNumberVector.size()+1, -1, RunNumberVector.size()));

	mOutList->Add(new TH2F("h_QA_vertex_x_TOF", "Vertex_x_position_TOFvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_y_TOF", "Vertex_y_position_TOFvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size())); 
	mOutList->Add(new TH2F("h_QA_vertex_z_TOF", "Vertex_z_position_TOFvsRunIndex", 100,-0.5, 0.5, RunNumberVector.size()+1, -1, RunNumberVector.size()));

//___BEMC QA histograms_______________________________________________________________________________________________________
	mOutList->Add(new TH2F("h_QA_BEMC_TOWId", "BEMC_Mached_tower_id_vs_RunIndex", 4801, 0, 4800, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_BSMD_nEta", "BSMD_Eta_wires_vs_RunIndex", 11, 0, 10, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_QA_BSMD_nPhi", "BSMD_Phi_wires_vs_RunIndex", 11, 0, 10, RunNumberVector.size()+1, -1, RunNumberVector.size()));
	mOutList->Add(new TH2F("h_pOverE", "E/p_vs_RunIndex", 200, 0, 10, RunNumberVector.size()+1, -1, RunNumberVector.size()));



//___dEdx_dNdx_nSigma_QA______________________________________________________________________________________________________





//-------------------------------------
  
	//histoInit(mOutFileBaseName, true); //for createQA()

//	set RefMultCorr
	mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt"); //for new StRefMultCorr, Run16, SL16j


   // -------------- USER VARIABLES -------------------------

	

	
//---Set branches------------------------------------------------------------------------------
	


//-----------------------Set branches END------------------------------------------------------------------------------

	mRunNumber = 0;
  return kStOK;
}

// _________________________________________________________
void StPicoDpmAnaMaker::ClearHF(Option_t *opt="") {
  return;
}

// _________________________________________________________
int StPicoDpmAnaMaker::FinishHF() {
   if( isMakerMode() != StPicoHFMaker::kWrite )
    //ntp_Dmeson->Write();
   //closeFile(); // for QA (analysis)
  return kStOK;
}

// _________________________________________________________
int StPicoDpmAnaMaker::MakeHF() {
  // -- process event
  //    ADD YOUR PROCESSING CODE HERE
  //    ... it is usefull to use the methods below
  //     - createCandidates()
  //     - analyzeCandidates()
//	cout<<"start"<<endl;
  std::clock_t start1 = std::clock();//kvapil
  if (isMakerMode() == StPicoHFMaker::kWrite) {
    createCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kRead) {
    // -- the reading back of the perviously written trees happens in the background
    analyzeCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kAnalyze) {
    //createCandidates(); //do not need fot QA of Run16
    //analyzeCandidates();
    //createQA();
  }

	RunId = mPicoDst->event()->runId();
	int RunIndex = -1; //default value for RunIndex (does not correspond to any RunId)

	for(unsigned int i=0; i<RunNumberVector.size();i++) //find corresponding RunIndex to a given RunId
	{
		if(RunNumberVector.at(i)==RunId)
		{
			RunIndex = i;
			break; //do not continue if coresponding RunId is found
		}
	}


   TH2F *h_mh1Cent = static_cast<TH2F*>(mOutList->FindObject("h_mh1Cent"));
   TH2F *h_mh1CentWg = static_cast<TH2F*>(mOutList->FindObject("h_mh1CentWg"));
   TH2F *h_mh1gRefmultCor = static_cast<TH2F*>(mOutList->FindObject("h_mh1gRefmultCor"));
   TH2F *h_mh1gRefmultCorWg = static_cast<TH2F*>(mOutList->FindObject("h_mh1gRefmultCorWg"));
   TH3F *h_mh2CentVz = static_cast<TH3F*>(mOutList->FindObject("h_mh2CentVz"));
   TH3F *h_mh2CentVzWg = static_cast<TH3F*>(mOutList->FindObject("h_mh2CentVzWg"));

		TH2D *h_QA_Vz = static_cast<TH2D*>(mOutList->FindObject("h_QA_Vz"));
		TH2D *h_QA_VzmVzVPD = static_cast<TH2D*>(mOutList->FindObject("h_QA_VzmVzVPD"));

		TH2D *h_QA_ZDC_rate = static_cast<TH2D*>(mOutList->FindObject("h_QA_ZDC_rate"));
		TH2D *h_QA_BBC_rate = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate"));

		TH2D *h_QA_ZDC_rate_pileUp = static_cast<TH2D*>(mOutList->FindObject("h_QA_ZDC_rate_pileUp"));
		TH2D *h_QA_BBC_rate_pileUp = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pileUp"));

		TH2D *h_QA_ZDC_rate_pileUp_TOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_ZDC_rate_pileUp_TOF"));
		TH2D *h_QA_BBC_rate_pileUp_TOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_BBC_rate_pileUp_TOF"));

		TH2F *h_QA_reweight_isNaN = static_cast<TH2F*>(mOutList->FindObject("h_QA_reweight_isNaN"));
		TH1I *h_QA_nEvents = static_cast<TH1I*>(mOutList->FindObject("h_QA_nEvents"));

		TH2D *h_QA_nTracks = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks"));
		TH2D *h_QA_nTracks_TPC = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_TPC"));
		TH2D *h_QA_nTracks_HFT = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT"));
		TH2D *h_QA_nTracks_HFT_PXL1 = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_PXL1"));
		TH2D *h_QA_nTracks_HFT_PXL2 = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_PXL2"));
		TH2D *h_QA_nTracks_HFT_IST = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_IST"));
		TH2D *h_QA_nTracks_HFT_SSD = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_SSD"));
		TH2D *h_QA_nTracks_HFT_IST_or_SSD = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_IST_or_SSD"));
		TH2D *h_QA_nTracks_TOF = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_TOF"));
		TH2D *h_QA_nTracks_BEMC = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_BEMC"));

		TH2D *h_QA_nTracks_HFT_etaCut = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut"));
		TH2D *h_QA_nTracks_HFT_etaCut_PXL1 = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_PXL1"));
		TH2D *h_QA_nTracks_HFT_etaCut_PXL2 = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_PXL2"));
		TH2D *h_QA_nTracks_HFT_etaCut_IST = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_IST"));
		TH2D *h_QA_nTracks_HFT_etaCut_SSD = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_SSD"));
		TH2D *h_QA_nTracks_HFT_etaCut_IST_or_SSD = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_HFT_etaCut_IST_or_SSD"));

	//Number of tracks vs. pT - get from integral of pT spectrum

//____general QA histograms (all tracks within (TPC, track quality) cuts, NO TOF and HFT)__________________________________
		TH2F *h_QA_pT = static_cast<TH2F*>(mOutList->FindObject("h_QA_pT"));
		TH2F *h_QA_eta = static_cast<TH2F*>(mOutList->FindObject("h_QA_eta"));
		TH2F *h_QA_phi = static_cast<TH2F*>(mOutList->FindObject("h_QA_phi"));

		TH2F *h_QA_vertex_x = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_x"));
		TH2F *h_QA_vertex_y = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_y"));
		TH2F *h_QA_vertex_z = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_z"));

		TH2F *h_QA_DCA_xy_TPC = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_xy_TPC"));
		TH2F *h_QA_DCA_xy_zoom_TPC = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_xy_zoom_TPC"));
		TH2F *h_QA_DCA_z_TPC = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_z_TPC"));
		TH2F *h_QA_DCA_z_zoom_TPC = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_z_zoom_TPC"));

//___HFT QA hitograms______________________________________________________________________________________________________
	
		TH2F *h_QA_pT_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_pT_HFT"));
		TH2F *h_QA_eta_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_eta_HFT"));
		TH2F *h_QA_phi_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_phi_HFT"));

		TH2F *h_QA_vertex_x_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_x_HFT"));
		TH2F *h_QA_vertex_y_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_y_HFT"));
		TH2F *h_QA_vertex_z_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_z_HFT"));

		TH2F *h_QA_DCA_xy_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_xy_HFT"));
		TH2F *h_QA_DCA_xy_zoom_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_xy_zoom_HFT"));
		TH2F *h_QA_DCA_z_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_z_HFT"));
		TH2F *h_QA_DCA_z_zoom_HFT = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_z_zoom_HFT"));

//___HFT QA hitograms, eta cut______________________________________________________________________________________________________
		TH2D *h_QA_nTracks_TPC_etaCut = static_cast<TH2D*>(mOutList->FindObject("h_QA_nTracks_TPC_etaCut"));

		TH2F *h_QA_pT_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_pT_HFT_etaCut"));
		TH2F *h_QA_eta_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_eta_HFT_etaCut"));
		TH2F *h_QA_phi_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_phi_HFT_etaCut"));

		TH2F *h_QA_vertex_x_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_x_HFT_etaCut"));
		TH2F *h_QA_vertex_y_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_y_HFT_etaCut"));
		TH2F *h_QA_vertex_z_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_z_HFT_etaCut"));

		TH2F *h_QA_DCA_xy_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_xy_HFT_etaCut"));
		TH2F *h_QA_DCA_xy_zoom_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_xy_zoom_HFT_etaCut"));
		TH2F *h_QA_DCA_z_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_z_HFT_etaCut"));
		TH2F *h_QA_DCA_z_zoom_HFT_etaCut = static_cast<TH2F*>(mOutList->FindObject("h_QA_DCA_z_zoom_HFT_etaCut"));


//____TOF QA histograms______________________________________________________________________________________________________
		TH2F *h_QA_pT_TOF = static_cast<TH2F*>(mOutList->FindObject("h_QA_pT_TOF"));
		TH2F *h_QA_eta_TOF = static_cast<TH2F*>(mOutList->FindObject("h_QA_eta_TOF"));
		TH2F *h_QA_phi_TOF = static_cast<TH2F*>(mOutList->FindObject("h_QA_phi_TOF"));

		TH2F *h_QA_vertex_x_TOF = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_x_TOF"));
		TH2F *h_QA_vertex_y_TOF = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_y_TOF"));
		TH2F *h_QA_vertex_z_TOF = static_cast<TH2F*>(mOutList->FindObject("h_QA_vertex_z_TOF"));

//___BEMC QA histograms_______________________________________________________________________________________________________
	TH2F *h_QA_BEMC_TOWId = static_cast<TH2F*>(mOutList->FindObject("h_QA_BEMC_TOWId"));
	TH2F *h_QA_BSMD_nEta = static_cast<TH2F*>(mOutList->FindObject("h_QA_BSMD_nEta"));
	TH2F *h_QA_BSMD_nPhi = static_cast<TH2F*>(mOutList->FindObject("h_QA_BSMD_nPhi"));
	TH2F *h_pOverE = static_cast<TH2F*>(mOutList->FindObject("h_pOverE"));


    StThreeVectorF pVtx = mPicoDst->event()->primaryVertex();
/*
//------------------SET TITLES OF HISTOGRAMS--------------------------------------------------------------------

//___detector, centrality and statistics (No. of events, tracks...) histograms
	h_mh1Cent->SetTitle(Form("EventsVsCentrality;cent;Counts_%d", RunId));
  h_mh1CentWg->SetTitle(Form("EventsVsCentrality;cent;Counts_%d", RunId));
  h_mh1gRefmultCor->SetTitle(Form("gRefmultCor;gRefmult;Counts_%d", RunId));
  h_mh1gRefmultCorWg->SetTitle(Form("gRefmultCorWg;gRefmultCorWg;Counts_%d", RunId));
  h_mh2CentVz->SetTitle(Form("CentralityVsVz;cent;Vz_%d", RunId));
  h_mh2CentVzWg->SetTitle(Form("CentralityVsVzWg;cent;Vz_%d", RunId));

	h_QA_reweight_isNaN->SetTitle(Form("reweight_isNaN_%d", RunId)); //check for bad weights for refmutCorr
	h_QA_nEvents->SetTitle(Form("Number_of_events_%d", RunId)); //number of events in run

	h_QA_nTracks->SetTitle(Form("Number_of_tracks_%d", RunId));	//total nuber of tracks in event (no cuts)
	h_QA_nTracks_TPC->SetTitle(Form("Number_of_TPC_tracks_%d", RunId)); //number of TPC tracks in run (heve to pass TPC cuts)
	h_QA_nTracks_HFT->SetTitle(Form("Number_of_HFT_tracks_%d", RunId)); //number of HFT tracks (have to pass track->isHFTTrack())
	h_QA_nTracks_HFT_PXL1->SetTitle(Form("Number_of_HFT_PXL1_tracks_%d", RunId)); //number of PXL1 tracks (have to pass track->hasPxl1Hit())
	h_QA_nTracks_HFT_PXL2->SetTitle(Form("Number_of_HFT_PXL2_tracks_%d", RunId)); //number of PXL2 tracks (have to pass track->hasPxl2Hit())
	h_QA_nTracks_HFT_IST->SetTitle(Form("Numbh_QA_ZDC_rateer_of_HFT_IST_tracks_%d", RunId));	//number of IST tracks (have to pass track->hasIstHit())
	h_QA_nTracks_HFT_SSD->SetTitle(Form("Number_of_HFT_SSD_tracks_%d", RunId));	//number of SSD tracks (have to pass track->hasSstHit())
	h_QA_nTracks_HFT_IST_or_SSD->SetTitle(Form("Number_of_HFT_IST_or_SSD_tracks_%d", RunId)); //number of IST or SST tracks (have to pass track->hasIstHit() || track->hasSstHit())
	h_QA_nTracks_TOF->SetTitle(Form("Number_of_TOF_tracks_%d", RunId)); //nuber of TOF tracks (have to have TOF info)

	//Number of tracks vs. pT - get from integral of pT spectrum

//____general QA histograms (all tracks within (TPC, track quality) cuts, NO TOF and HFT)__________________________________
	h_QA_pT->SetTitle(Form("Transverse_momentum_TPC_%d", RunId)); //check binning and range
	h_QA_eta->SetTitle(Form("Pseudorapidity_TPC_%d", RunId));
	h_QA_phi->SetTitle(Form("Azimuthal_angle_TPC_%d", RunId));

	h_QA_vertex_x->SetTitle(Form("Vertex_x_position_%d", RunId)); 
	h_QA_vertex_y->SetTitle(Form("Vertex_y_position_%d", RunId)); 
	h_QA_vertex_z->SetTitle(Form("Vertex_z_position_%d", RunId));

	h_QA_DCA_xy_TPC->SetTitle(Form("DCA_xy_TPC_%d", RunId));
	h_QA_DCA_xy_zoom_TPC->SetTitle(Form("DCA_xy_zoom_TPC_%d", RunId));
	h_QA_DCA_z_TPC->SetTitle(Form("DCA_z_TPC_%d", RunId));
	h_QA_DCA_z_zoom_TPC->SetTitle(Form("DCA_z_zoom_TPC_%d", RunId));


//___HFT QA hitograms______________________________________________________________________________________________________
	h_QA_pT_HFT->SetTitle(Form("Transverse_momentum_HFT_%d", RunId)); //check binning and range
	h_QA_eta_HFT->SetTitle(Form("Pseudorapidity_HFT", RunId));
	h_QA_phi_HFT->SetTitle(Form("Azimuthal_angle_HFT_%d", RunId));

	h_QA_vertex_x_HFT->SetTitle(Form("Vertex_x_position_HFT_%d", RunId)); 
	h_QA_vertex_y_HFT->SetTitle(Form("Vertex_h_QA_ZDC_ratey_position_HFT_%d", RunId)); 
	h_QA_vertex_z_HFT->SetTitle(Form("Vertex_z_position_HFT_%d", RunId));

	h_QA_DCA_xy_HFT->SetTitle(Form("DCA_xy_HFT_%d", RunId));
	h_QA_DCA_xy_zoom_HFT->SetTitle(Form("DCA_xy_zoom_HFT_%d", RunId));
	h_QA_DCA_z_HFT->SetTitle(Form("DCA_z_HFT_%d", RunId));
	h_QA_DCA_z_zoom_HFT->SetTitle(Form("DCA_z_zoom_HFT_%d", RunId));

//____TOF QA histograms______________________________________________________________________________________________________
	h_QA_pT_TOF->SetTitle(Form("Transverse_momentum_TOF_%d", RunId)); //check binning and range
	h_QA_eta_TOF->SetTitle(Form("Pseudorapidity_TOF_%d", RunId));
	h_QA_phi_TOF->SetTitle(Form("Azimuthal_angle_TOF_%d", RunId));

	h_QA_vertex_x_TOF->SetTitle(Form("Vertex_x_position_TOF_%d", RunId)); 
	h_QA_vertex_y_TOF->SetTitle(Form("Vertex_y_position_TOF_%d", RunId)); 
	h_QA_vertex_z_TOF->SetTitle(Form("Vertex_z_position_TOF_%d", RunId));
*/

    mRefmultCorrUtil->init(mPicoDst->event()->runId());

      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;

      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;
      int const centrality = mRefmultCorrUtil->getCentralityBin9();
      const double reweight = mRefmultCorrUtil->getWeight();
      const double refmultCor = mRefmultCorrUtil->getRefMultCorr();

			const int zero = 0;
			const int one = 1;

			if(std::isnan(reweight))
			{
					h_QA_reweight_isNaN->Fill(one, RunIndex);					
			}
			else
			{
					h_QA_reweight_isNaN->Fill(zero, RunIndex); //0 - reweight OK, 1 - not OK (reweight is NaN)
			}

      h_mh1gRefmultCor->Fill(refmultCor, RunIndex);
      h_mh1gRefmultCorWg->Fill(refmultCor, RunIndex, reweight);
      h_mh1Cent->Fill(centrality, RunIndex);
      h_mh1CentWg->Fill(centrality, RunIndex, reweight);
      h_mh2CentVz->Fill(centrality, pVtx.z(), RunIndex);
      h_mh2CentVzWg->Fill(centrality, pVtx.z(), RunIndex, reweight);		

			h_QA_nEvents->Fill(RunIndex); //number of 0 filled to this histogram = number of events	

			h_QA_ZDC_rate->Fill(mPicoDst->event()->ZDCx()/1000., RunIndex);
			h_QA_BBC_rate->Fill(mPicoDst->event()->BBCx()/1000., RunIndex);		
		
			float vertex_x_QA = pVtx.x();
			float vertex_y_QA = pVtx.y();
			float vertex_z_QA = pVtx.z();

			h_QA_Vz->Fill(vertex_z_QA, RunIndex);
			h_QA_VzmVzVPD->Fill(fabs(vertex_z_QA - mPicoDst->event()->vzVpd()), RunIndex);

    UInt_t nTracks = mPicoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
        StPicoTrack const* trk = mPicoDst->track(iTrack);
        if (!trk) continue;

				StPhysicalHelixD helix = trk->helix(mPicoDst->event()->bField()); //SL16j, Vanek
        StThreeVectorF momentum = trk->gMom(pVtx, mPicoDst->event()->bField());
				float dca_xy_QA = float(helix.geometricSignedDistance(pVtx.x(), pVtx.y())); //dca_xy
				StThreeVectorF dcaPoint = helix.at(helix.pathLength(vertex_x_QA, vertex_y_QA));
        float dca_z_QA = dcaPoint.z() - vertex_z_QA; //check

				float pT_QA = trk->gPt();
				float eta_QA = momentum.pseudoRapidity();
				float phi_QA = momentum.phi();

				float Beta = mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()); //SL16j, Vanek
	      if(!isnan(Beta) && Beta > 0){
					h_QA_ZDC_rate_pileUp_TOF->Fill(mPicoDst->event()->ZDCx()/1000., RunIndex);
					h_QA_BBC_rate_pileUp_TOF->Fill(mPicoDst->event()->BBCx()/1000., RunIndex);
				}			

				h_QA_nTracks->Fill(pT_QA, RunIndex); 

				//add all relevant TPC and track quality cuts here				
				if(!(mHFCuts->hasGoodNHitsFitMinHist(trk))) continue;
				if(!(mHFCuts->hasGoodNHitsFitnHitsMax(trk))) continue;
        if(!(mHFCuts->hasGoodEta(momentum))) continue;
				if(!(mHFCuts->hasGoodPtQA(trk))) continue;
				if(!(mHFCuts->hasGoodDCAMaxGlob(helix, pVtx))) continue;

				h_QA_ZDC_rate_pileUp->Fill(mPicoDst->event()->ZDCx()/1000., RunIndex);
				h_QA_BBC_rate_pileUp->Fill(mPicoDst->event()->BBCx()/1000., RunIndex);
				
				//Fill all TPC histograms here
				
				h_QA_nTracks_TPC->Fill(pT_QA, RunIndex); 
				
				h_QA_pT->Fill(pT_QA, RunIndex);
				h_QA_eta->Fill(eta_QA, RunIndex);
				h_QA_phi->Fill(phi_QA, RunIndex);
	
				h_QA_vertex_x->Fill(vertex_x_QA, RunIndex);
				h_QA_vertex_y->Fill(vertex_y_QA, RunIndex);
				h_QA_vertex_z->Fill(vertex_z_QA, RunIndex);
		
				h_QA_DCA_xy_TPC->Fill(dca_xy_QA, RunIndex);
				h_QA_DCA_xy_zoom_TPC->Fill(dca_xy_QA, RunIndex);
				h_QA_DCA_z_TPC->Fill(dca_z_QA, RunIndex);
				h_QA_DCA_z_zoom_TPC->Fill(dca_z_QA, RunIndex);

				//---------------FILL HFT QA INFORMATION-------------------------------------------------------
				if(trk->isHFTTrack())
				{
					//Fill all HFT histograms here	

					h_QA_nTracks_HFT->Fill(pT_QA, RunIndex);

					h_QA_pT_HFT->Fill(pT_QA, RunIndex);
					h_QA_eta_HFT->Fill(eta_QA, RunIndex);
					h_QA_phi_HFT->Fill(phi_QA, RunIndex);

					h_QA_vertex_x_HFT->Fill(vertex_x_QA, RunIndex);
					h_QA_vertex_y_HFT->Fill(vertex_y_QA, RunIndex);
					h_QA_vertex_z_HFT->Fill(vertex_z_QA, RunIndex);
		
					h_QA_DCA_xy_HFT->Fill(dca_xy_QA, RunIndex);
					h_QA_DCA_xy_zoom_HFT->Fill(dca_xy_QA, RunIndex);
					h_QA_DCA_z_HFT->Fill(dca_z_QA, RunIndex);
					h_QA_DCA_z_zoom_HFT->Fill(dca_z_QA, RunIndex);

				}

				if(trk->hasPxl1Hit())
				{
					//Fill all PXL1 histograms here

					h_QA_nTracks_HFT_PXL1->Fill(pT_QA, RunIndex);

				}

				if(trk->hasPxl2Hit())
				{
					//Fill all PXL2 histograms here

					h_QA_nTracks_HFT_PXL2->Fill(pT_QA, RunIndex);

				}

				if(trk->hasIstHit())
				{
 					//Fill all IST histograms here

					h_QA_nTracks_HFT_IST->Fill(pT_QA, RunIndex);

				}

				if(trk->hasSstHit())
 				{
					//Fill all SSD histograms here

					h_QA_nTracks_HFT_SSD->Fill(pT_QA, RunIndex);

				} 

				if( (trk->hasIstHit() || trk->hasSstHit() ))
				{ 
					//Fill all (IST or SSD) histograms here

					h_QA_nTracks_HFT_IST_or_SSD->Fill(pT_QA, RunIndex);

				}
				
				 //---------------FILL HFT QA INFORMATION, DIFFERENT ETA CUT-- -----------------------------------------------------
				if(mHFCuts->hasGoodEtaHFT(momentum))
				{
					h_QA_nTracks_TPC_etaCut->Fill(pT_QA, RunIndex); //for HFT_etaCut histograms
				}

				if( trk->isHFTTrack() && mHFCuts->hasGoodEtaHFT(momentum))
	 			{
					//Fill all HFT histograms here
					h_QA_nTracks_HFT_etaCut->Fill(pT_QA, RunIndex); 

					h_QA_pT_HFT_etaCut-> Fill(pT_QA, RunIndex);
					h_QA_eta_HFT_etaCut->Fill(eta_QA, RunIndex);
					h_QA_phi_HFT_etaCut->Fill(phi_QA, RunIndex);

					h_QA_vertex_x_HFT_etaCut->Fill(vertex_x_QA, RunIndex);
					h_QA_vertex_y_HFT_etaCut->Fill(vertex_y_QA, RunIndex);
					h_QA_vertex_z_HFT_etaCut->Fill(vertex_z_QA, RunIndex);
		
					h_QA_DCA_xy_HFT_etaCut->Fill(dca_xy_QA, RunIndex);
					h_QA_DCA_xy_zoom_HFT_etaCut->Fill(dca_xy_QA, RunIndex);
					h_QA_DCA_z_HFT_etaCut->Fill(dca_z_QA, RunIndex);
					h_QA_DCA_z_zoom_HFT_etaCut->Fill(dca_z_QA, RunIndex);

				}

				if(trk->hasPxl1Hit() && mHFCuts->hasGoodEtaHFT(momentum))
				{
					//Fill all PXL1 histograms here

					h_QA_nTracks_HFT_etaCut_PXL1->Fill(pT_QA, RunIndex);

				}

				if(trk->hasPxl2Hit() && mHFCuts->hasGoodEtaHFT(momentum))
				{
					//Fill all PXL2 histograms here

					h_QA_nTracks_HFT_etaCut_PXL2->Fill(pT_QA, RunIndex);

				}

				if(trk->hasIstHit() && mHFCuts->hasGoodEtaHFT(momentum))
				{
					//Fill all IST histograms here

					h_QA_nTracks_HFT_etaCut_IST->Fill(pT_QA, RunIndex);

				}

				if(trk->hasSstHit() && mHFCuts->hasGoodEtaHFT(momentum))
				{
					//Fill all SSD histograms here

					h_QA_nTracks_HFT_etaCut_SSD->Fill(pT_QA, RunIndex);

				}

				if( (trk->hasIstHit() || trk->hasSstHit() ) && mHFCuts->hasGoodEtaHFT(momentum))
				{
					//Fill all (IST or SSD) histograms here

					h_QA_nTracks_HFT_etaCut_IST_or_SSD->Fill(pT_QA, RunIndex);

				}

 				//------------------------------FILL TOF QA INFORMATION---------------------------------------       
				float piBeta = mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()); //SL16j, Vanek
        Int_t piTofAvailable = 1; //piTofAvailable = 0 -> OK; = 1 -> no TOF info
        if(!isnan(piBeta) && piBeta > 0){
					piTofAvailable = 0;
					//Fill all TOF histograms here

					h_QA_nTracks_TOF->Fill(pT_QA, RunIndex);

					h_QA_pT_TOF->Fill(pT_QA, RunIndex);
					h_QA_eta_TOF->Fill(eta_QA, RunIndex);
					h_QA_phi_TOF->Fill(phi_QA, RunIndex);
	
					h_QA_vertex_x_TOF->Fill(vertex_x_QA, RunIndex);
					h_QA_vertex_y_TOF->Fill(vertex_y_QA, RunIndex);
					h_QA_vertex_z_TOF->Fill(vertex_z_QA, RunIndex);

   			}	
				

				//-------------------------------------FILL BEMC QA INFORMATION------------------------

				if(trk->bemcPidTraitsIndex()>-1) //check good BEMC tracks - no BEMC info i current SL16j Run16 PicoDst production
				{
					StPicoBEmcPidTraits * Emc =  mPicoDst->bemcPidTraits(trk->bemcPidTraitsIndex()); //check that this works		
					//cout<<"BEMC track"<<endl;
					h_QA_nTracks_BEMC->Fill(pT_QA, RunIndex);		
				
					h_QA_BEMC_TOWId->Fill(Emc->btowId(), RunIndex);
					h_QA_BSMD_nEta->Fill(Emc->bemcSmdNEta(), RunIndex);
					h_QA_BSMD_nPhi->Fill(Emc->bemcSmdNPhi(), RunIndex);
					h_pOverE->Fill(momentum.mag()/Emc->bemcE0(), RunIndex);
				}
				

				
					
 

   
  } // .. end tracks loop


 
  return kStOK;
}

int StPicoDpmAnaMaker::createQA(){
       //int const currentRun = mPicoHFEvent->runId();
       //if(currentRun != mRunNumber)
     //  {
       //mRunNumber = currentRun;
      mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;
//      cout<<"Q1"<<endl;
//      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
//      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
//	  cout<<"Q2"<<endl;
      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;

       int const centrality = mRefmultCorrUtil->getCentralityBin9();
       const double reweight = mRefmultCorrUtil->getWeight();
       const double refmultCor = mRefmultCorrUtil->getRefMultCorr();
       //mHists->addCent(refmultCor, centrality, reweight, pVtx.z());
       UInt_t nTracks = mPicoDst->numberOfTracks();


       for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
       {
          StPicoTrack const* trk = mPicoDst->track(iTrack);
          if (!trk) continue;
          //StPhysicalHelixD helix = trk->helix(); //SL16d
		  		StPhysicalHelixD helix = trk->helix(mPicoDst->event()->bField()); //SL16j, Vanek
          float dca = float(helix.geometricSignedDistance(mPrimVtx));
          StThreeVectorF momentum = trk->gMom(mPrimVtx, mPicoDst->event()->bField());

					StThreeVectorF dcaPoint = helix.at(helix.pathLength(mPrimVtx.x(), mPrimVtx.y()));
          float dcaZ = dcaPoint.z() - mPrimVtx.z();
          double dcaXy = helix.geometricSignedDistance(mPrimVtx.x(), mPrimVtx.y());

					


					// Saving to TTree
					// save all variables (branches) defined in InitHF()

	  

					//ntp_Dmeson->Fill()



/*
         	// if (!isGoodQaTrack(trk, momentum, dca)) continue; pt, nhits, pseudorap
					if (!(mHFCuts->hasGoodPtQA(trk))) continue;
					if (!(mHFCuts->hasGoodNHitsFitMinHist(trk))) continue;
        	if (!(mHFCuts->hasGoodEta(momentum))) continue;

          StThreeVectorF dcaPoint = helix.at(helix.pathLength(mPrimVtx.x(), mPrimVtx.y()));
          float dcaZ = dcaPoint.z() - mPrimVtx.z();
          double dcaXy = helix.geometricSignedDistance(mPrimVtx.x(), mPrimVtx.y());

          bool tpcPion = false;
          bool tpcKaon = false;
          bool tpcProton = false;
		  		if(mHFCuts->hasGoodTPCnSigmaPion(trk)) tpcPion = true;
		  		if(mHFCuts->hasGoodTPCnSigmaKaon(trk)) tpcKaon = true;
	      	if(mHFCuts->hasGoodTPCnSigmaProton(trk)) tpcProton = true;
          //float hBeta = mHFCuts->getTofBetaBase(trk); //SL16d
		  		float hBeta = mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()); //SL16j, Vanek
          bool hTofAvailable = !isnan(hBeta) && hBeta > 0;

          bool tofPion = false;
          bool tofKaon = false;
          bool tofProton = false;

	      	if(fabs(1./hBeta - sqrt(1+M_PION_PLUS*M_PION_PLUS/(momentum.mag()*momentum.mag())))<=mHFCuts->getCutTOFDeltaOneOverBeta(StHFCuts::kPion)) tofPion = true;
          if(fabs(1./hBeta - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(momentum.mag()*momentum.mag())))<=mHFCuts->getCutTOFDeltaOneOverBeta(StHFCuts::kKaon)) tofKaon = true;
		  		if(fabs(1./hBeta - sqrt(1+M_PROTON*M_PROTON/(momentum.mag()*momentum.mag())))) tofProton = true;

          bool goodPion = (hTofAvailable && tofPion && tpcPion) || (!hTofAvailable && tpcPion);//Always require TPC
          bool goodKaon = (hTofAvailable && tofKaon && tpcKaon) || (!hTofAvailable && tpcKaon);
          bool goodProton = (hTofAvailable && tofProton && tpcProton) || (!hTofAvailable && tpcProton);

          if (trk  && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton)){ //createQA() not used in this case - see myDpmAnalysisQA
             addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //add Dca distribution
          }
          if (trk  && fabs(dca) < 1.5 && (goodPion || goodKaon || goodProton)){
             //std::cout<<"1: "<<goodPion<<" "<< goodKaon<<" "<<  goodProton<<" "<<  momentum.perp()<<" "<<  centrality<<" "<<  momentum.pseudoRapidity()<<" "<<  momentum.phi()<<" "<<  mPrimVtx.z()<<std::endl;
             addTpcDenom1(goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add Tpc Denominator
          }
          if (trk && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton) && fabs(dcaXy) < 1. && fabs(dcaZ) < 1.){
             addHFTNumer1(goodPion, goodKaon, goodProton, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add HFT Numerator
          }
*/
       } // .. end tracks loop
   return 0;
}

// _________________________________________________________
int StPicoDpmAnaMaker::createCandidates() {
  // Creating candidates for D+- 3 body decay
  // D+- -> K+2Pi decay

  for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
    StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
    // -- Pion selection      

    for (unsigned short idxPion2 = idxPion1+1; idxPion2 < mIdxPicoPions.size(); ++idxPion2) {
      StPicoTrack const *pion2 = mPicoDst->track(mIdxPicoPions[idxPion2]);
      // -- Pion selection
      if ( !isCloseTracks(pion1,pion2,mPrimVtx, mBField)) continue; 

      for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
        StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
        // -- Kaon selection
        // -- TOF
		//if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon) ) continue; //SL16d
        if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon, mPicoDst->event()->bField()), StHFCuts::kKaon) ) continue; //SL16j, Vanek
        if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1]|| mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion2] || mIdxPicoPions[idxPion1] == mIdxPicoPions[idxPion2]) continue;
				if ( !isCloseTracks(pion1,kaon,mPrimVtx, mBField)) continue; 
				if ( !isCloseTracks(kaon,pion2,mPrimVtx, mBField)) continue; 
        // -- Making triplet
        StHFTriplet triplet(pion1,pion2,kaon,mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mIdxPicoPions[idxPion1],mIdxPicoPions[idxPion2],mIdxPicoKaons[idxKaon], mPrimVtx, mBField);
        if(mHFCuts->hasGoodTripletdV0Max(triplet)) continue;
        if (!mHFCuts->isGoodSecondaryVertexTriplet(triplet)) continue;
        mPicoHFEvent->addHFSecondaryVertexTriplet(&triplet);

      }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion2 = idxPion1+1; idxPion2 < mIdxPicoPions.size(); ++idxPion2)
  } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)
 return kStOK;
}

// _________________________________________________________
int StPicoDpmAnaMaker::analyzeCandidates() {

//not used in Run16 QA

 return kStOK;
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
  // -- good hadron
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isPion(StPicoTrack const * const trk) const {
  // -- good pion
   StThreeVectorF t = trk->pMom();
   if (fabs(t.pseudoRapidity()) > 1.) return false; //pridano fabs 1212
   //if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion) ) return false; //SL16d
	 if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kPion) ) return false; //SL16j, Vanek
   if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isKaon(StPicoTrack const * const trk) const {
  // -- good kaon
  StThreeVectorF t = trk->pMom();
  if (fabs(t.pseudoRapidity()) > 1.) return false;//pridano fabs 1212
	//if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false; //SL16d
  if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kKaon) ) return false; //SL16j, Vanek
  if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
} 

// _________________________________________________________
bool StPicoDpmAnaMaker::isProton(StPicoTrack const * const trk) const {
  // -- good proton
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}

double StPicoDpmAnaMaker::DCA(StPicoTrack const * const trk, StThreeVectorF const & vtx) const {
  // -- particle DCA
/*  StPhysicalHelixD pHelix = trk->dcaGeometry().helix(); //SL16d
  pHelix.moveOrigin(pHelix.pathLength(vtx)); 
  return ((pHelix.origin() - vtx).mag());
*/
	return ((trk->origin() - vtx).mag()); //SL16j, Vanek
}


bool StPicoDpmAnaMaker::isCloseTracks(StPicoTrack const * const trk1, StPicoTrack const * const trk2, StThreeVectorF const & vtx, float bField) const {
/* SL16d  
  StPhysicalHelixD p1Helix = trk1->dcaGeometry().helix();
  StPhysicalHelixD p2Helix = trk2->dcaGeometry().helix();
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  if( ( p1Helix.origin()-vtx ).mag()>0.2 || ( p2Helix.origin()-vtx ).mag()>0.2 ) return false;
*/
	if( ( trk1->origin()-vtx ).mag()>0.2 || ( trk2->origin()-vtx ).mag()>0.2 ) return false; //SL16j, Vanek

  //Requires loading constants
//  StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss); //SL16d
//  StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);

	StThreeVectorF const p1Mom = trk1->gMom(); //SL16j, Vanek
  StThreeVectorF const p2Mom = trk2->gMom();
  StPhysicalHelixD const p1StraightLine(p1Mom, trk1->origin(), 0, trk1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, trk2->origin(), 0, trk2->charge());
  //DCA
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);
  float const dca = (p1AtDcaToP2-p2AtDcaToP1).mag();
  if(dca > 0.009) return false;
// -- good pair
  return true;
}

//-----------------------------------------------------------------------------

void StPicoDpmAnaMaker::histoInit(TString fileBaseName, bool fillQaHists){
  TString m_ParticleName[m_nParticles] = {"Pion", "Kaon", "Proton"};

   float m_EtaEdgeDca[m_nEtasDca+1] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0}; //replace bottom!!!
   float m_PhiEdgeDca[m_nPhisDca + 1] = {-3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};
   float m_VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};//replace bottom!!!
   float m_CentEdgeDca[m_nCentsDca + 1] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
   float m_PtEdgeDca[m_nPtsDca + 1] = {0.3, 0.4, 0.5, 0.6,  0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 6.0, 12.0};
   float m_EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0}; //replace bottom!!!
   float m_PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};//replace bottom!!!
   float m_VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};//replace bottom!!!
   float m_CentEdgeRatio[m_nCentsRatio + 1] = { -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
   float m_PtEdgeRatio[m_nPtsRatio + 1] =
   {
      0.3, 0.4, 0.5, 0.6 , 0.7 , 0.8 , 0.9 ,
      1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
      2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 ,
      3.4 , 3.8 , 4.2 , 4.6 , 5.0 ,  5.5 ,
      6. , 6.5 , 7.0 , 8.0 , 9.0 , 10. , 11,  12.0
   };
  float m_DcaEdgeDca[m_nDcasDca + 1] =
   {
     -1 , -0.96 , -0.92 , -0.88 , -0.84 , -0.8 , -0.76 , -0.72 , -0.68 , -0.64 , -0.6 , -0.56 , -0.52 , -0.48 , -0.44 , -0.4 , -0.36 , -0.32 , -0.28 , -0.24 , -0.2 , -0.16 , -0.12 ,  -0.08,
     -0.078 , -0.075 , -0.072 , -0.069 , -0.066 , -0.063 , -0.06 , -0.057 , -0.054 , -0.051 , -0.048 , -0.045 , -0.042 , -0.039 , -0.036 , -0.033 , -0.03 , -0.027 , -0.024 , -0.021 , -0.018 , -0.015 , -0.012 ,
      -0.01 , -0.0096 , -0.0092 , -0.0088 , -0.0084 , -0.008 , -0.0076 , -0.0072 , -0.0068 , -0.0064 , -0.006 , -0.0056 , -0.0052 , -0.0048 , -0.0044 , -0.004 , -0.0036 , -0.0032 , -0.0028 , -0.0024 , -0.002 , -0.0016 , -0.0012 , -0.0008 , -0.0004 , 0 , 0.0004 , 0.0008 , 0.0012 , 0.0016 , 0.002 , 0.0024 , 0.0028 , 0.0032 , 0.0036 , 0.004 , 0.0044 , 0.0048 , 0.0052 , 0.0056 , 0.006 , 0.0064 , 0.0068 , 0.0072 , 0.0076 , 0.008 , 0.0084 , 0.0088 , 0.0092 , 0.0096 , 0.01 ,
      0.012 , 0.015 , 0.018 , 0.021 , 0.024 , 0.027 , 0.03 , 0.033 , 0.036 , 0.039 , 0.042 , 0.045 , 0.048 , 0.051 , 0.054 , 0.057 , 0.06 , 0.063 , 0.066 , 0.069 , 0.072 , 0.075 , 0.078 ,
      0.08 , 0.12 , 0.16 , 0.2 , 0.24 , 0.28 , 0.32 , 0.36 , 0.4 , 0.44 , 0.48 , 0.52 , 0.56 , 0.6 , 0.64 , 0.68 , 0.72 , 0.76 , 0.8 , 0.84 , 0.88 , 0.92 , 0.96 , 1
   };



   //set in private
   // for(int temp = 0;temp<m_nParticles;temp++) m_ParticleName[temp]=temp_ParticleName[temp];
	//for(int temp = 0;temp<m_nEtasDca+1;temp++) m_EtaEdgeDca[temp]=temp_EtaEdgeDca[temp];
	//for(int temp2 = 0;temp2<m_nPhisDca+1;temp2++) m_PhiEdgeDca[temp2]=temp_PhiEdgeDca[temp2];
	/*for(int temp = 0;temp<m_nVzsDca+1;temp++) m_VzEdgeDca[temp]=temp_VzEdgeDca[temp];
	for(int temp = 0;temp<m_nCentsDca+1;temp++) m_CentEdgeDca[temp]=temp_CentEdgeDca[temp];
	for(int temp = 0;temp<m_nPtsDca+1;temp++) m_PtEdgeDca[temp]=temp_PtEdgeDca[temp];
	for(int temp = 0;temp<m_nEtasRatio+1;temp++) m_EtaEdgeRatio[temp]=temp_EtaEdgeRatio[temp];
	for(int temp = 0;temp<m_nPhisRatio+1;temp++) m_PhiEdgeRatio[temp]=temp_PhiEdgeRatio[temp];
	for(int temp = 0;temp<m_nVzsRatio+1;temp++) m_VzEdgeRatio[temp]=temp_VzEdgeRatio[temp];
	for(int temp = 0;temp<m_nCentsRatio+1;temp++) m_CentEdgeRatio[temp]=temp_CentEdgeRatio[temp];
	for(int temp = 0;temp<m_nPtsRatio+1;temp++) m_PtEdgeRatio[temp]=temp_PtEdgeRatio[temp];
	for(int temp = 0;temp<m_nDcasDca+1;temp++) m_DcaEdgeDca[temp]=temp_DcaEdgeDca[temp];*/

   mFillQaHists = fillQaHists;
   mOutFile = new TFile(fileBaseName+".hists.root", "RECREATE");
   //mOutFile = new TFile(Form("%s.hists.root", fileBaseName.Data()), "RECREATE");
/*
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasDca; iEta++){
         for (int iVz = 0; iVz < m_nVzsDca; iVz++){
            for (int iCent = 0; iCent < m_nCentsDca; iCent++){
               mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent] = NULL;
            }
         }
      }
   }

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++){
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++){
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++){
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
            }
         }
      }
   }*/

  

   TH1::SetDefaultSumw2();
   if (!mFillQaHists) return;
/*
   mh1Cent         = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5); //histograms already initialized and filled in InitHF() and MakeHF()??
   mh1CentWg         = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1gRefmultCor  = new TH1F("mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700);
   mh1gRefmultCorWg  = new TH1F("mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700);
   mh2CentVz         = new TH2F("mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh2CentVzWg = new TH2F("mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
*/
   //Add some HFT ratio plots
   mh2Tpc1PtCent  = new TH2F("mh2Tpc1PtCent", "Tpc tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2HFT1PtCent  = new TH2F("mh2HFT1PtCent", "HFT tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2Tpc1PhiVz  = new TH2F("mh2Tpc1PhiVz", "Tpc tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm
   mh2HFT1PhiVz  = new TH2F("mh2HFT1PhiVz", "HFT tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++){
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++){
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++){
             
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2Tpc1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), "mh2Tpc1PtCent_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", m_EtaEdgeRatio[iEta], m_VzEdgeRatio[iVz], m_PhiEdgeRatio[iPhi]), m_nPtsRatio, m_PtEdgeRatio, m_nCentsRatio, m_CentEdgeRatio); //Dca 1.cm
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2HFT1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), "mh2HFT1PtCent_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", m_EtaEdgeRatio[iEta], m_VzEdgeRatio[iVz], m_PhiEdgeRatio[iPhi]), m_nPtsRatio, m_PtEdgeRatio, m_nCentsRatio, m_CentEdgeRatio); //Dca 1.cm
            }
         }
      }
   }

   // Add some Dca, resolution
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasDca; iEta++){
         for (int iVz = 0; iVz < m_nVzsDca; iVz++){
            for (int iCent = 0; iCent < m_nCentsDca; iCent++){
               
   	    	   mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent]  = new TH3F(Form("mh3DcaXyZPtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iCent),"mh3DcaXyZPt_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Cent%2.1f;p_{T}(GeV/c);DcaXy(cm);DcaZ(cm)", m_EtaEdgeDca[iEta], m_VzEdgeDca[iVz], m_CentEdgeDca[iCent]), m_nPtsDca, m_PtEdgeDca, m_nDcasDca, m_DcaEdgeDca, m_nDcasDca, m_DcaEdgeDca); //Dca 1.cm
            }
         }
      }
   }

   mh3DcaPtCent  = new TH3F("mh3DcaPtCent", "mh3DcaPtCent;p_{T}(GeV/c);cent;Dca(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaXyPtCent  = new TH3F("mh3DcaXyPtCent", "mh3DcaXyPtCent;p_{T}(GeV/c);cent;DcaXy(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaZPtCent  = new TH3F("mh3DcaZPtCent", "mh3DcaZPtCent;p_{T}(GeV/c);cent;DcaZ(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm

}

//-----------------------------------------------------------------------
void StPicoDpmAnaMaker::addTpcDenom1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
      if(EtaIndex == -1) return;
   if(PhiIndex == -1) return;
   if(VzIndex == -1) return;
   //std::cout<<"2: "<<IsPion<<" "<<IsKaon<<" "<<IsProton<<" "<<pt<<" "<<centrality<<" "<<Eta<<" "<<Phi<<" "<<Vz<<" "<<EtaIndex<<" "<<PhiIndex<<" "<<VzIndex<<std::endl;
   
   if (IsPion){
      mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      //if(mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]) std::cout<<"true"<<<<std::endl;
      //std::cout<<pt<<" "<<centrality<<std::endl;
   }
   if (IsKaon){
      mh2Tpc1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsProton){
      mh2Tpc1PtCentPartEtaVzPhi[2][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2Tpc1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2Tpc1PhiVz->Fill(Phi, Vz);
}
//-----------------------------------------------------------------------
void StPicoDpmAnaMaker::addHFTNumer1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
   if(EtaIndex == -1) return;
   if(PhiIndex == -1) return;
   if(VzIndex == -1) return;
   if (IsPion){
      mh2HFT1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsKaon){
      mh2HFT1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsProton){
      mh2HFT1PtCentPartEtaVzPhi[2][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2HFT1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2HFT1PhiVz->Fill(Phi, Vz);
}
//---------------------------------------------------------------------
void StPicoDpmAnaMaker::addDcaPtCent(float dca, float dcaXy, float dcaZ, bool IsPion, bool IsKaon, bool IsProton, float pt,  int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexDca(Eta);
   int VzIndex = getVzIndexDca(Vz);
   if(EtaIndex == -1) return;
   if(VzIndex == -1) return;

   if (centrality < 0) return; // remove bad centrality, only keep 9 centralities
   if (IsPion){
      mh3DcaXyZPtCentPartEtaVzPhi[0][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsKaon){
      mh3DcaXyZPtCentPartEtaVzPhi[1][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsProton){
      mh3DcaXyZPtCentPartEtaVzPhi[2][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   mh3DcaPtCent->Fill(pt, centrality, dca);
   mh3DcaXyPtCent->Fill(pt, centrality, dcaXy);
   mh3DcaZPtCent->Fill(pt, centrality, dcaZ);
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getEtaIndexDca(float Eta){
   float EtaEdgeDca[m_nEtasDca+1] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0};
   for (int i = 0; i < m_nEtasDca; i++){
	 if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
         return i;
   }
   //std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nEtasDca -1;
   return -1;
}

//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getVzIndexDca(float Vz){
  float VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};
   for (int i = 0; i < m_nVzsDca; i++){
      if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nVzsDca - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getEtaIndexRatio(float Eta){
  float EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0}; 
   for (int i = 0; i < m_nEtasRatio; i++){
      if ((Eta >= EtaEdgeRatio[i]) && (Eta < EtaEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nEtasRatio - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getPhiIndexRatio(float Phi){
  float PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};
   for (int i = 0; i < m_nPhisRatio; i++){
      if ((Phi >= PhiEdgeRatio[i]) && (Phi < PhiEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
  // return m_nPhisRatio - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getVzIndexRatio(float Vz){
  float VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};
   for (int i = 0; i < m_nVzsRatio; i++) {
      if ((Vz >= VzEdgeRatio[i]) && (Vz < VzEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBLE WRONG"<<std::endl;
  // return m_nVzsRatio - 1;
  return -1;
}

void StPicoDpmAnaMaker::addCent(const double refmultCor, int centrality, const double reweight, const float vz)
{
   mh1gRefmultCor->Fill(refmultCor);
   mh1gRefmultCorWg->Fill(refmultCor, reweight);
   mh1Cent->Fill(centrality);
   mh1CentWg->Fill(centrality, reweight);
   mh2CentVz->Fill(centrality, vz);
   mh2CentVzWg->Fill(centrality, vz, reweight);
}

//---------------------------------------------------------------------
void StPicoDpmAnaMaker::closeFile()
{
   mOutFile->cd();

   mh1Cent->Write();
   mh1CentWg->Write();
   mh1gRefmultCor->Write();
   mh1gRefmultCorWg->Write();
   mh2CentVz->Write();
   mh2CentVzWg->Write();

   //HFT ratio QA
   mh2Tpc1PtCent->Write();
   mh2Tpc1PhiVz->Write();
   mh2HFT1PhiVz->Write();
   mh2HFT1PtCent->Write();

   //HFT DCA Ratio
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < m_nEtasDca; iEta++)
      {
         for (int iVz = 0; iVz < m_nVzsDca; iVz++)
         {
            for (int iCent = 0; iCent < m_nCentsDca; iCent++)
            {
               mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent]->Write();
            }
         }
      }
   }
  // std::cout<<"tuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu"<<m_nParticles<<" "<<m_nEtasRatio<<std::endl;

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++)
      {
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++)
         {
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
            }
         }
      }
   }

   mh3DcaPtCent->Write();
   mh3DcaXyPtCent->Write();
   mh3DcaZPtCent->Write();

   // nt->Write();
   mOutFile->Write();
   mOutFile->Close();
   //mOutFile->Delete();
}



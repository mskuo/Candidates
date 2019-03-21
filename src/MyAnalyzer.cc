// This is a file for select four-vectors from particles in datasets
//
//
// system include files
#include <memory>

// user include files, general
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//extra header files 
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include <string>

//for root histogramming
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

//for TNtuples 
#include "TNtuple.h"

//for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

//for PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//for muon collection
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//for gsfelectron collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

//for vertex information
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//for PFJet collection
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//for HLTrigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

//class declaration

class DemoAnalyzer: public edm::EDAnalyzer 
{
	public:
		explicit DemoAnalyzer(const edm::ParameterSet&);
		~DemoAnalyzer();

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		bool providesGoodLumisection( const edm::Event& iEvent) ;
	TFile *myfile ;	
	//declare ttree
	TTree *mytree ;

	//declare variables
	int 	eventno,	nparticles,	type[3000],	njets,		nvertex,	charge[3000] ;
	double 	px[3000],	py[3000],	pz[3000]	;
	double 	pt[3000],	eta[3000],	phi[3000],	e[3000] ;
	double  j_pt[3000],	j_eta[3000],	j_phi[3000],	j_e[3000] ;
	double  x[50],		y[50],		z[50] ;

	//variable for HLTrigger
//	#define N_TRIGGER_BOOKINGS 5788
	int trgBook[9] ;
	int trgCount ;	
	int n_trg ;

	const std::string TriggerBooking[9] = {
	"HLT_Mu13_Mu8_v17",	//000
	"HLT_Mu17_Mu8_v17",	//001
	"HLT_Mu5_v18",		//002
	"HLT_Mu8_v16",		//003
	"HLT_Mu12_v16",		//004
	"HLT_Mu17_v3",		//005
	"HLT_Mu24_v14",		//006
	"HLT_Mu30_v14",		//007
	"HLT_Mu40_v12"		//008
	} ;

};


DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig) 
{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs ;

	//intialize the varibles
	eventno = 9999 ; nparticles = 9999 ; nvertex = 9999 ; njets = 9999 ; 

}


DemoAnalyzer::~DemoAnalyzer() 
{
	// do anything here that needs to be done at destruction time
	// (eg. close files, deallocate resources etc.)
}


//method called for each event
void DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
	using namespace edm ;
	using namespace reco ;
	using namespace std ;

	#ifdef THIS_IS_AN_EVENT_EXAMPLE
		Handle<ExampleData> pIn ;
		iEvent.getByLabel ( "example" , pIn) ;
	#endif

	#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
		ESHandle<SetupData> pSetup ;
		iSetup.get<SetupRecord>().get(pSetup) ;
	#endif

	//Event is to be analyzed
	
	LogInfo("Demo") 
	<< "Starting to analyze \n"
	<< "Event number: " << (iEvent.id()).event()
	<< "\nRun number: " << iEvent.run()
	<< "\nLumisection: " << iEvent.luminosityBlock() ;

	Handle<reco::PFCandidateCollection> particles ;
	iEvent.getByLabel( "particleFlow" , particles ) ;
	
	Handle<reco::VertexCollection> primvtxHandle ;
	iEvent.getByLabel("offlinePrimaryVertices" , primvtxHandle ) ;

	Handle<reco::PFJetCollection>  jets ;
	iEvent.getByLabel( "ak5PFJets" , jets ) ;


	//HLT
	edm::Handle<TriggerResults> TrgResultsHandle ; //catch triggerresults
	const pat::helper::TriggerMatchHelper matchHelper ;

	edm::InputTag trigResultsTag( "TriggerResults" , "" , "HLT" ) ;
	bool with_TriggerResults = iEvent.getByLabel( trigResultsTag , TrgResultsHandle ) ;
	if ( !with_TriggerResults )
	{
		std::cout << "Sorry there is no TriggerResult in the file" << std::endl ; 
	}
	else
	{
		//get the names of the triggers
		const edm::TriggerNames &TrgNames = iEvent.triggerNames( *TrgResultsHandle ) ;
		trgCount = 0 ;
		n_trg = TrgNames.size() ;

		//Check the trigger list in the event
//		for ( unsigned k = 0 ; k < 440 ; k++ )
//		{
//			cout << TrgNames.triggerNames()[k] << endl ;
//		}
		for ( int i = 0 ; i < 9 ; i++ )
		{
			unsigned int TrgIndex = TrgNames.triggerIndex( TriggerBooking[i] ) ;
//			cout << "i = " <<  i << endl ;
//			cout << "TriggerBooking[i] = " << TriggerBooking[i] << endl ;
//			cout << "TrgNames.size() = " << TrgNames.size()  << endl ;
//			cout << "TrgIndex 	 = " << TrgIndex << endl ;

			if ( TrgIndex == TrgNames.size() )
			{
				trgBook[i] = -4 ; // The trigger path is not known in this event.
			}
			else if ( !TrgResultsHandle -> wasrun( TrgIndex ) )
			{
				trgBook[i] = -3 ; // The trigger path was not included in this event.
			}
			else if ( !TrgResultsHandle -> accept( TrgIndex ) )
			{
				trgBook[i] = -2 ; // The trigger path was not accepted in this event.
			}
			else if (  TrgResultsHandle -> error ( TrgIndex ) )
			{
				trgBook[i] = -1 ; // The trigger path has an error in this event.
			}
			else
			{
				trgBook[i] = +1 ; // It's triggered.
				trgCount++ ;
			}
		}
		//save all trigger
//		EvtInfo.nHLT = TrgNames.size() ;
//		for ( unsigned int i = 0 ; i < TrgNames.size() ; i++ )
//		{
//			EvtInfo.hltBits[i] = ( TrgResultsHandle -> accept(i) == true ) ? 1:0 ;
//		}
	}//end(!with_TriggerResults)







	//initialze the counting variables
	int c = 0 ;
	int jc = 0 ;

	eventno = (iEvent.id()).event() ;
	nparticles = (particles -> size())-1 ; 
	njets = (jets -> size())-1 ;
	
	//vertex
	reco::VertexCollection primvtx ;

	if ( primvtxHandle.isValid() )
	{
		primvtx = *primvtxHandle ;
	}
	else
	{
		edm::LogInfo("Demo") << "No primary vertex available from EventSetup \n" ;
	}

	nvertex = primvtx.size() ;

	for ( int i = 0 ; i != nvertex ; i++ )
	{
		math::XYZPoint point(primvtx[i].position()) ;
		x[i] = point.x() ;
		y[i] = point.y() ;
		z[i] = point.z() ;
	}		

	//Loop over all the particles of current Event
	for (reco::PFCandidateCollection::const_iterator it = particles -> begin() ; it != particles -> end() ; it++ )
	{ 			
		type[c] = it -> particleId() ;
		px[c]  = it -> px() ; 
		py[c]  = it -> py() ;
		pz[c]  = it -> pz() ;
		pt[c]  = it -> pt() ;
		eta[c] = it -> eta() ;
		phi[c] = it -> phi() ;
		e[c]  = it -> energy() ;
		charge[c] = it -> charge() ;

		c++ ;
	}
	
	for (reco::PFJetCollection::const_iterator it = jets -> begin() ; it != jets -> end() ; it++ )
	{ 			
		j_pt[jc]  = it -> pt() ;
		j_eta[jc] = it -> eta() ;
		j_phi[jc] = it -> phi() ;
		j_e[jc]  = it -> energy() ;

		jc++ ;
	}


	mytree -> Fill() ;	

}

void DemoAnalyzer::beginJob()
{
	myfile = new TFile("myfile.root" , "recreate");
	//TTree
	mytree = new TTree ( "Candidates" , "Candidates" ) ;
	mytree -> Branch ( "eventno" , &eventno , "eventno/I" ) ;
	mytree -> Branch ( "nparticles" , &nparticles , "nparticles/I" ) ;
	mytree -> Branch ( "type" , type , "type[nparticles]/I" ) ;
	mytree -> Branch ( "px" , px , "px[nparticles]/D" ) ;
	mytree -> Branch ( "py" , py , "py[nparticles]/D" ) ;
	mytree -> Branch ( "pz" , pz , "pz[nparticles]/D" ) ;
	mytree -> Branch ( "pt" , pt , "pt[nparticles]/D" ) ;
	mytree -> Branch ( "eta" , eta , "eta[nparticles]/D" ) ;
	mytree -> Branch ( "phi" , phi , "phi[nparticles]/D" ) ;
	mytree -> Branch ( "e" , e , "e[nparticles]/D" ) ;
	mytree -> Branch ( "charge" , charge , "charge[nparticles]/I" ) ;
	
	mytree -> Branch ( "nvertex" , &nvertex , "nvertex/I" ) ;
	mytree -> Branch ( "x" , x , "x[nvertex]/D" ) ;
	mytree -> Branch ( "y" , y , "y[nvertex]/D" ) ;
	mytree -> Branch ( "z" , z , "z[nvertex]/D" ) ; 

	mytree -> Branch ( "njets" , &njets , "njets/I" ) ;
	mytree -> Branch ( "j_pt" , j_pt , "j_pt[njets]/D" ) ;
	mytree -> Branch ( "j_eta" , j_eta , "j_eta[njets]/D" ) ;
	mytree -> Branch ( "j_phi" , j_phi , "j_phi[njets]/D" ) ;
	mytree -> Branch ( "j_e" , j_e , "j_e[njets]/D" ) ;

	mytree -> Branch ( "trgBook" , trgBook , "trgBook[9]/I" ) ;
	mytree -> Branch ( "n_trg" , &n_trg , " n_trg/I" ) ;
}


void DemoAnalyzer::endJob()
{
	mytree -> Write() ;
	myfile -> Close() ;
}

DEFINE_FWK_MODULE(DemoAnalyzer);	

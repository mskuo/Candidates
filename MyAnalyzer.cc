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
	int eventno, nparticles, type[3000] ;
//	double px[3000], py[3000], pz[3000] ;
	double pt[3000], eta[3000], phi[3000] ;
	int nvertex ;
	double x[50], y[50], z[50] ;
};


DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig) 
{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs ;

	//intialize the varibles
	eventno = 9999 ; nparticles = 9999 ; nvertex = 9999 ;
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

	//initialze the counting variables
	int c = 0 ;

	eventno = (iEvent.id()).event() ;
	nparticles = (particles -> size())-1 ; 

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
//		px[c]  = it -> px() ;
//		py[c]  = it -> py() ;
//		pz[c]  = it -> pz() ;
		pt[c]  = it -> pt() ;
		eta[c] = it -> eta() ;
		phi[c] = it -> phi() ;
		c++ ;
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
//	mytree -> Branch ( "px" , px , "px[nparticles]/D" ) ;
//	mytree -> Branch ( "py" , py , "py[nparticles]/D" ) ;
//	mytree -> Branch ( "pz" , pz , "pz[nparticles]/D" ) ;
	mytree -> Branch ( "pt" , pt , "pt[nparticles]/D" ) ;
	mytree -> Branch ( "eta" , eta , "eta[nparticles]/D" ) ;
	mytree -> Branch ( "phi" , phi , "phi[nparticles]/D" ) ;
	
	mytree -> Branch ( "nvertex" , &nvertex , "nvertex/I" ) ;
	mytree -> Branch ( "x" , x , "x[nvertex]/D" ) ;
	mytree -> Branch ( "y" , y , "y[nvertex]/D" ) ;
	mytree -> Branch ( "z" , z , "z[nvertex]/D" ) ; 
}

void DemoAnalyzer::endJob()
{
	mytree -> Write() ;
	myfile -> Close() ;
}

DEFINE_FWK_MODULE(DemoAnalyzer);	

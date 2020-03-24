// Package:    MuonAnalyzer
//              AOD version
// 
/**\class 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <vector>
#include "TTree.h"
#include <string>
#include <iostream>
#include "TMath.h"
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "L1Trigger/L1TNtuples/interface/MuonID.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
/*namespace edm {
  class ConfigurationDescriptions;
  }*/

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
template<typename T1>
class MuonAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;
  typedef edm::AssociationMap<edm::OneToValue<std::vector<T1>, float > > T1IsolationMap;

public:
  explicit MuonAODAnalyzer(const edm::ParameterSet&);
  ~MuonAODAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void genMu(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;
 

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken l1MuonsToken_;
  edm::EDGetToken tracksToken_;
  edm::EDGetToken GenToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;
  std::vector<std::string> HLTPaths;
  std::vector<std::string> HLTFilters;

 
 unsigned int event,run_number,ls,nmuons,ndimuons,ntrkmuons,ntracks;
 std::vector<float> pvertex_x,pvertex_y,pvertex_z;
 float beam_x=0,beam_y=0,beam_z=0;

//gen
 std::vector<float> genmuon_pt,genmuon_eta,genmuon_phi,genmuon_charge,genmuon_momId,genmuon_grandmomId,genmuon_grandgrandmomId;
//trg
 std::vector<bool> hltpaths;
 std::vector<float> l1muon_pt,l1muon_eta,l1muon_phi,hltmuon_pt,hltmuon_eta,hltmuon_phi;
 std::vector<unsigned int> hltmuon_pathId;
 unsigned int hlt_firedMuindex;
 float  hlt_firedMuDr;
 //muon
 std::vector<bool> muon_soft, muon_medium, muon_loose, muon_tight, 
                   muon_global, muon_pf, muon_vtxcand;
 std::vector<float> muon_pt, muon_eta, muon_phi, muon_charge, muon_vx, muon_vy,
                    muon_vz, muon_iso, muon_dxy, muon_edxy, muon_dz, muon_edz,
                    muon_normChi2, muon_nvalidhits, muon_nmatchstations,
                    muon_npixelhits, muon_trklayers;
 //jpsi
 //
 std::vector<float> dimuon_pt,dimuon_eta,dimuon_phi,dimuon_ctxy,dimuon_x,dimuon_y,dimuon_z,dimuon_prob,dimuon_mass,dimuon_mass_unfit,dimuon_mu1pt,dimuon_mu1eta,dimuon_mu1phi,dimuon_mu2pt,dimuon_mu2eta,dimuon_mu2phi,dimuon_cos2D,dimuon_lxy,dimuon_elxy;
 std::vector<unsigned int> dimuon_mu1id,dimuon_mu2id;
 //tracks
  std::vector<float> track_pt,track_phi,track_eta,track_charge,track_vx,track_vy,track_vz,track_dxy,track_edxy,track_dz,track_edz,track_trklayers,track_pxllayers,track_normchi2;
  std::vector<bool> track_vtxcand;
 //trkmu
 std::vector<float> trkmuon_pt,trkmuon_eta,trkmuon_phi,trkmuon_ctxy,trkmuon_x,trkmuon_y,trkmuon_z,trkmuon_prob,trkmuon_mass,trkmuon_mass_unfit,trkmuon_mupt,trkmuon_mueta,trkmuon_muphi,trkmuon_trkpt,trkmuon_trketa,trkmuon_trkphi,trkmuon_cos2D,trkmuon_lxy,trkmuon_elxy;
 std::vector<unsigned int> trkmuon_muid,trkmuon_trkid;
 //  edm::ParameterSet const& conf;
 

  edm::Service<TFileService> fs;
  TTree * t1;
  bool MuMuVtx=false; bool MuTrkVtx=false; bool SaveTracks=true;
  bool data=true;
  bool SkipEvtWithNoFire=false; 
  bool CutSelection=false;
  bool TrgSelection=false;
  bool GenSelection=true;
  double TrgMuMatchCone=100;
  //muon sel
  double minMuPtCut=0;
  double maxMuTrgDzCut=100; 
  double maxMuTrgDrCut=100;
  bool SoftQCut=false;
  //track sel
  double TrkPtCut=-1;
  double DzTrkTrg=-1;
  //gen sel
  double GenRecoMatchCone=100;
  double MomId=0;
  //vtx cuts - before fit
  double minDiMuPtCut=0; 
  double McutMin=0; 
  double McutMax=100;
  double minMuTrkPtCut=0;
  double DzMuMu=-1; 
  double DzMuTrk=-1;
  //vtx cuts - after fit
  bool SkipEvtNoHLTmuMatch=false; double maxDrMuMatch=-1;
  bool SkipEvtNoMuVtx=false; bool SkipEvtNoTrkVtx=false;
  bool SaveMatchGenTrkOnly=false;
  double ProbCut=-1; double CosCut=-1.1;
     // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
template<typename T1>
MuonAODAnalyzer<T1>::MuonAODAnalyzer(const edm::ParameterSet& iConfig): 
 beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"))),
 vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons"))),
  tracksToken_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  GenToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen"))),
trgresultsToken_(consumes<edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("triggerresults"))),
 trigobjectsToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
 HLTPaths(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),
 HLTFilters(iConfig.getParameter<std::vector<std::string>>("HLTFilters"))
{
  edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
 data=runParameters.getParameter<bool>("IsData");
 SkipEvtWithNoFire=runParameters.getParameter<bool>("SkipEvtWithNoFire");
 SkipEvtNoMuVtx=runParameters.getParameter<bool>("SkipEvtNoMuVtx");
 SkipEvtNoTrkVtx=runParameters.getParameter<bool>("SkipEvtNoTrkVtx");
 SkipEvtNoHLTmuMatch=runParameters.getParameter<bool>("SkipEvtNoHLTmuMatch");
 SaveTracks=runParameters.getParameter<bool>("SaveTracks");
 SaveMatchGenTrkOnly=runParameters.getParameter<bool>("SaveMatchGenTrkOnly");
 CutSelection=runParameters.getParameter<bool> ("CutSelection");
 TrgSelection=runParameters.getParameter<bool> ("TrgSelection");
 GenSelection=runParameters.getParameter<bool> ("GenSelection");

 TrgMuMatchCone=runParameters.getParameter<double>("TrgMuMatchCone");
 minMuPtCut=runParameters.getParameter<double>("MuPtCut");
 maxMuTrgDzCut=runParameters.getParameter<double>("DzMuTrgCut");
 maxMuTrgDrCut=runParameters.getParameter<double>("DrMuTrgCut");
 SoftQCut=runParameters.getParameter<bool>("SoftQCut");
 
 minDiMuPtCut=runParameters.getParameter<double>("minDiMuPtCut");
 McutMin=runParameters.getParameter<double>("MresonanceMin");
 McutMax=runParameters.getParameter<double>("MresonanceMax");
 DzMuTrk=runParameters.getParameter<double>("DzMuTrk");
 DzMuMu=runParameters.getParameter<double>("DzMuMu");
 DzTrkTrg=runParameters.getParameter<double>("DzTrkTrg");
 TrkPtCut=runParameters.getParameter<double>("TrkPtCut");
 minMuTrkPtCut=runParameters.getParameter<double>("minMuTrkPtCut");
 
 GenRecoMatchCone=runParameters.getParameter<double>("GenRecoMatchCone");
 MomId=runParameters.getParameter<double>("MomId");
 MuMuVtx=runParameters.getParameter<bool>("MuMuVtx");
 MuTrkVtx=runParameters.getParameter<bool>("MuTrkVtx");
 
 ProbCut=runParameters.getParameter<double>("MuTrkProbCut");
 CosCut=runParameters.getParameter<double>("MuTrkCosCut");
}

template<typename T1>
MuonAODAnalyzer<T1>::~MuonAODAnalyzer()
{
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
   // do anything here that needs to be done at desctruction time

}


//
// member functions
//

// ------------ method called for each event  ------------

template<typename T1>
void MuonAODAnalyzer<T1>::genMu(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<std::vector<reco::GenParticle>> genPart;
  iEvent.getByToken(GenToken_,genPart); 
  for ( const reco::GenParticle  & gen: *genPart){
    if (fabs(gen.pdgId())!=13) continue;
    if (gen.status()!=1) continue;
    if (gen.numberOfMothers()==0) continue;
    genmuon_pt.push_back(gen.pt()); genmuon_eta.push_back(gen.eta());
    genmuon_phi.push_back(gen.phi()); genmuon_charge.push_back(gen.charge());  
    const reco::Candidate * mom=gen.mother();
    genmuon_momId.push_back(mom->pdgId());
    if (mom->numberOfMothers()>0)
      genmuon_grandmomId.push_back(mom->mother()->pdgId());
    else 
      genmuon_grandmomId.push_back(-999);
    if (mom->numberOfMothers()>0 && mom->mother()->numberOfMothers()>0)
      genmuon_grandgrandmomId.push_back(mom->mother()->mother()->pdgId());
    else 
      genmuon_grandgrandmomId.push_back(999);
  }
}


  
template<typename T1>
void
MuonAODAnalyzer<T1>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  //clean vectors
  //gen
  genmuon_pt.clear(),genmuon_eta.clear(),genmuon_phi.clear(),genmuon_charge.clear(),genmuon_momId.clear(),genmuon_grandmomId.clear(); genmuon_grandgrandmomId.clear();
 //trg
 l1muon_pt.clear(),l1muon_eta.clear(),l1muon_phi.clear();
 hltmuon_pt.clear(),hltmuon_eta.clear(),hltmuon_phi.clear(),hltmuon_pathId.clear();
 hltpaths.clear();
 //reco muon
 muon_soft.clear(),muon_medium.clear(),muon_loose.clear(),muon_tight.clear();
 muon_global.clear(); muon_pf.clear(); muon_vtxcand.clear();
 muon_pt.clear(),muon_eta.clear(),muon_phi.clear(),muon_charge.clear();
 muon_vx.clear(),muon_vy.clear(),muon_vz.clear(),muon_iso.clear();
 muon_dxy.clear(); muon_edxy.clear(); muon_dz.clear(); muon_edz.clear();
 muon_normChi2.clear(); muon_nvalidhits.clear(); muon_nmatchstations.clear();
 muon_npixelhits.clear(); muon_trklayers.clear();
 //track
 track_pt.clear(); track_phi.clear(); track_eta.clear(); track_charge.clear(); 
 track_vx.clear(); track_vy.clear(); track_vz.clear(); track_dxy.clear(); 
 track_edxy.clear(); track_dz.clear(); track_edz.clear(); track_trklayers.clear();
 track_pxllayers.clear(); track_normchi2.clear(); track_vtxcand.clear();
 //dimu
 dimuon_pt.clear(),dimuon_eta.clear(),dimuon_phi.clear(),dimuon_ctxy.clear();
 dimuon_x.clear(),dimuon_y.clear(),dimuon_z.clear(),dimuon_prob.clear(),dimuon_mass.clear();
 dimuon_mass_unfit.clear(),dimuon_mu1pt.clear(),dimuon_mu1eta.clear(),dimuon_mu1phi.clear();
 dimuon_mu2pt.clear(),dimuon_mu2eta.clear(),dimuon_mu2phi.clear();
 dimuon_cos2D.clear(); dimuon_lxy.clear(); dimuon_elxy.clear();
 dimuon_mu1id.clear(); dimuon_mu2id.clear();
//trkmu
 trkmuon_pt.clear(),trkmuon_eta.clear(),trkmuon_phi.clear(),trkmuon_ctxy.clear();
 trkmuon_x.clear(),trkmuon_y.clear(),trkmuon_z.clear(),trkmuon_prob.clear(),trkmuon_mass.clear();
 trkmuon_mass_unfit.clear(),trkmuon_mupt.clear(),trkmuon_mueta.clear(),trkmuon_muphi.clear();
 trkmuon_trkpt.clear(),trkmuon_trketa.clear(),trkmuon_trkphi.clear();
 trkmuon_cos2D.clear(); trkmuon_lxy.clear(); trkmuon_elxy.clear();
 trkmuon_muid.clear(); trkmuon_trkid.clear();
 //vtxs
 pvertex_x.clear(); pvertex_y.clear(); pvertex_z.clear();
 //zero
 ndimuons=0; ntrkmuons=0; ntracks=0; nmuons=0;
 //Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot); 
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  //continue if there are no vertices
  if (vertices->size()==0) return;
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);
  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_,tracks);
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  edm::Handle<l1t::MuonBxCollection> l1Muons;
  iEvent.getByToken(l1MuonsToken_,l1Muons);
  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_ ,triggerObjects);
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::TriggerNames trigName;
  trigName = iEvent.triggerNames(*trigResults);

  event++;
  //run stuff
  run_number=iEvent.id().run(); ls=iEvent.luminosityBlock();  
  beam_x= theBeamSpot->x0(); beam_y= theBeamSpot->y0();
  beam_z= theBeamSpot->z0();
  reco::TrackBase::Point  vertex_point; 
  const  reco::Vertex firstGoodVertex=vertices->front();
  for (const reco::Vertex &vtx : *vertices) {
    if (vtx.isFake() || !vtx.isValid ()) continue;
    pvertex_x.push_back(vtx.x()); 
    pvertex_y.push_back(vtx.y()); 
    pvertex_z.push_back(vtx.z());
  }
  if (pvertex_x.size()==0) return;
  vertex_point.SetCoordinates(pvertex_x[0],pvertex_y[0],pvertex_z[0]);
  
  std::vector< std::pair<float,float> > GensEtaPhi;

  if (!data){
    genMu(iEvent, iSetup);
    if (GenSelection){
      for (unsigned int igmuon=0; igmuon<genmuon_pt.size(); ++igmuon){
	
        if ( fabs(genmuon_momId[igmuon]) != MomId ) continue;
	
	GensEtaPhi.emplace_back( std::make_pair(genmuon_eta[igmuon],genmuon_phi[igmuon]) );
      }
    }
  }
  for(typename std::vector< l1t::Muon >::const_iterator mu=l1Muons->begin(0); mu !=l1Muons->end(0); mu++){
    l1muon_pt.push_back(mu->et()); l1muon_eta.push_back(mu->eta());
    l1muon_phi.push_back(mu->phi());
  }
  if( trigResults.failedToGet() ){
    if (SkipEvtWithNoFire) return;
    for (unsigned int trg=0; trg<HLTPaths.size(); trg++)
       hltpaths.push_back(false);
  }
  int N_Triggers = trigResults->size();
  bool EvtFire=false;
  for (unsigned int itrg=0; itrg<HLTPaths.size(); itrg++){
   bool fire=false;
   for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
     TString TrigPath =trigName.triggerName(i_Trig);
     if (!trigResults->accept(i_Trig)) continue;
     if (TrigPath.Contains(HLTPaths[itrg])){
       fire=true; EvtFire=true;
     }
   }
   hltpaths.push_back(fire);
  }

  if (!EvtFire && SkipEvtWithNoFire) return;



  trigger::TriggerObjectCollection allTriggerObjects = triggerObjects->getObjects(); 
 for (unsigned int ihlt=0; ihlt<HLTFilters.size(); ++ihlt){
   size_t filterIndex = (*triggerObjects).filterIndex(InputTag(HLTFilters[ihlt],"","HLT"));
    if (filterIndex < (*triggerObjects).sizeFilters()) {
      const trigger::Keys &keys = (*triggerObjects).filterKeys(filterIndex);
      for (size_t j = 0; j < keys.size(); j++) {
        trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
        if (fabs(foundObject.id())!=13) continue;
        hltmuon_pt.push_back(foundObject.pt());
        hltmuon_eta.push_back(foundObject.eta());
        hltmuon_phi.push_back(foundObject.phi());
        hltmuon_pathId.push_back(ihlt); 
      }
    }
 }

  float etamax=0,phimax=0,maxhlt=0;
  for (unsigned int ihlt=0; ihlt<hltmuon_pt.size(); ihlt++){
    if (hltmuon_pt[ihlt]<maxhlt) continue;
    maxhlt=hltmuon_pt[ihlt]; etamax=hltmuon_eta[ihlt]; phimax=hltmuon_phi[ihlt];
  }

      
 float minDR=100; hlt_firedMuindex=999;  hlt_firedMuDr=-999;
 reco::Muon trgmu;
 for (const reco::Muon & mu: *muons){
   if (minDR<deltaR(mu.eta(),mu.phi(),etamax,phimax)) continue;
   minDR=deltaR(mu.eta(),mu.phi(),etamax,phimax);
   hlt_firedMuindex=&mu-&muons->at(0); 
   hlt_firedMuDr=minDR;
   trgmu=mu;
 }

 if (SkipEvtNoHLTmuMatch && minDR>TrgMuMatchCone && TrgMuMatchCone>0) return;

 //save muons
 std::vector<reco::Muon> dimucands;  std::vector<unsigned int> dimuindex;
 for (const reco::Muon &mu : *muons){    
    muon_pt.push_back(mu.pt());   muon_phi.push_back(mu.phi());
    muon_eta.push_back(mu.eta()); muon_charge.push_back(mu.charge());    
    muon_vx.push_back(mu.vx());   muon_vy.push_back(mu.vy());
    muon_vz.push_back(mu.vz());   muon_medium.push_back(muon::isMediumMuon(mu));
    muon_loose.push_back(muon::isLooseMuon(mu));
    muon_tight.push_back(muon::isTightMuon(mu,firstGoodVertex));
    muon_soft.push_back(muon::isSoftMuon(mu,firstGoodVertex));
    const MuonPFIsolation&  isol=mu.pfIsolationR04();
    muon_iso.push_back((isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt());  
    muon_dxy.push_back(mu.bestTrack()->dxy(vertex_point));
    muon_dz.push_back(mu.bestTrack()->dz(vertex_point));
    muon_edxy.push_back(mu.dxyError());
    muon_edz.push_back(mu.dzError());    
    muon_global.push_back(mu.isGlobalMuon());
    muon_pf.push_back(mu.isPFMuon());
    if (mu.isGlobalMuon()){
       muon_normChi2.push_back(mu.globalTrack()->normalizedChi2()); 
       muon_nvalidhits.push_back(mu.globalTrack()->hitPattern()
                                 .numberOfValidMuonHits() );
    } else{
       muon_normChi2.push_back(mu.bestTrack()->normalizedChi2()); 
       muon_nvalidhits.push_back(mu.bestTrack()->hitPattern()
                                 .numberOfValidMuonHits() );
    }
    muon_nmatchstations.push_back(mu.numberOfMatchedStations());
    if(mu.innerTrack().isNonnull() ){
       muon_npixelhits.push_back(mu.innerTrack()->hitPattern()
                                 .numberOfValidPixelHits());  
       muon_trklayers.push_back(mu.innerTrack()->hitPattern()
                                .trackerLayersWithMeasurement());
    } else {
       muon_npixelhits.push_back(-1);  muon_trklayers.push_back(-1);
    }
    //fill candidates
    //first check gen mathed cands will run on mc only the cvector will be empty if we do not cfg it properly (possible only in mc)
    if (GenSelection){
      for ( auto & genEtaPhi: GensEtaPhi){
        if ( deltaR( genEtaPhi.first, genEtaPhi.second,
                    mu.eta(), mu.phi() ) > GenRecoMatchCone )
            continue;
        dimucands.emplace_back(mu); 
        dimuindex.push_back(nmuons);
        break;
      }
    }
    //check if there is an hlt muon and if use only hlt vertex selected
    if (TrgSelection && hlt_firedMuindex==nmuons){
       dimucands.emplace_back(mu); dimuindex.push_back(nmuons);
    }
    // else check if muon passes cuts
    if (CutSelection){
      if (mu.pt()>minMuPtCut &&  (fabs(mu.vz()-trgmu.vz())<maxMuTrgDzCut || maxMuTrgDzCut<0 ) &&  (deltaR(mu.eta(),mu.phi(),trgmu.eta(),trgmu.phi()) <maxMuTrgDrCut ||   maxMuTrgDrCut<0) && ( muon::isSoftMuon(mu,firstGoodVertex) || !SoftQCut)  ){
   	dimucands.emplace_back(mu);  dimuindex.push_back(nmuons);
     }
    }
    if ( dimuindex.size()>0 && nmuons == dimuindex.back()) 
       muon_vtxcand.push_back(true);
    else 
       muon_vtxcand.push_back(false);
    nmuons++;
  }


 std::vector<reco::Track> trkcands;  std::vector<unsigned int> trkindex;
 
 for (const reco::Track &trk : (*tracks)){    
   if (fabs(trk.eta())>2.5) continue;
   if (!trk.quality(Track::highPurity)) continue;
   if (CutSelection || TrgSelection){
     if(fabs(trk.vz()-trgmu.vz())>DzTrkTrg && DzTrkTrg>0) continue;
     if (trk.pt()<TrkPtCut) continue;
     trkcands.push_back(trk);  trkindex.push_back(ntracks);
    }  
    if (GenSelection){
      for ( auto & genEtaPhi: GensEtaPhi){
        if ( deltaR( genEtaPhi.first, genEtaPhi.second, 
                     trk.eta(), trk.phi() ) > GenRecoMatchCone ) 
            continue;
        trkcands.push_back(trk);  
        trkindex.push_back(ntracks);
        break;
     }
   }
    // save only if we want in ntuple
   if (SaveTracks || (SaveMatchGenTrkOnly && trkindex.size()>0 && ntracks==trkindex.back() ) ){
     track_pt.push_back(trk.pt()); track_phi.push_back(trk.phi());
     track_eta.push_back(trk.eta()); track_charge.push_back(trk.charge());    
     track_vx.push_back(trk.vx()); track_vy.push_back(trk.vy());
     track_vz.push_back(trk.vz()); track_dxy.push_back(trk.dxy(vertex_point));
     track_dz.push_back(trk.dz(vertex_point));
     track_edxy.push_back(trk.dxyError()); track_edz.push_back(trk.dzError()); 
     track_trklayers.push_back(trk.hitPattern().trackerLayersWithMeasurement());
     track_pxllayers.push_back(trk.hitPattern().pixelLayersWithMeasurement());
     track_normchi2.push_back(trk.normalizedChi2());  
    if ( trkindex.size()>0 && ntracks==trkindex.back() )
      track_vtxcand.push_back(true);
    else 
      track_vtxcand.push_back(false);
    ntracks++;
   }   
   
  }  

 
 if (MuMuVtx){  
  for (std::vector<reco::Muon>::iterator mu1=dimucands.begin(); mu1!=dimucands.end(); ++mu1){
   for (std::vector<reco::Muon>::iterator mu2=mu1+1; mu2!=dimucands.end(); ++mu2){
      if (mu1->charge()==mu2->charge()) continue;
      
      TLorentzVector vmu1,vmu2; 
      vmu1.SetPtEtaPhiM(mu1->pt(),mu1->eta(),mu1->phi(),0.105);
      vmu2.SetPtEtaPhiM(mu2->pt(),mu2->eta(),mu2->phi(),0.105);
      if ( (vmu1+vmu2).M()<McutMin || (vmu1+vmu2).M()>McutMax) continue;
      dimuon_mass_unfit.push_back((vmu1+vmu2).M());
      std::vector<reco::TransientTrack> dimutrk; dimutrk.reserve(2);
      dimutrk.emplace_back(reco::TransientTrack(*(mu1->bestTrack()),&(*bFieldHandle)));
      dimutrk.emplace_back(reco::TransientTrack(*(mu2->bestTrack()),&(*bFieldHandle)));
      KalmanVertexFitter vtxFitter(true);
      TransientVertex dimuvtx=vtxFitter.vertex(dimutrk);
      if(!dimuvtx.isValid()) continue; 
      dimuon_x.push_back(dimuvtx.position().x());
      dimuon_y.push_back(dimuvtx.position().y());
      dimuon_z.push_back(dimuvtx.position().z());
      dimuon_prob.push_back(ChiSquaredProbability(dimuvtx.totalChiSquared(), dimuvtx.degreesOfFreedom()));
      std::vector<reco::TransientTrack> refited=dimuvtx.refittedTracks();
      vmu1.SetPtEtaPhiM(refited[0].track().pt(),refited[0].track().eta(),refited[0].track().phi(),0.105);
      vmu2.SetPtEtaPhiM(refited[1].track().pt(),refited[1].track().eta(),refited[1].track().phi(),0.105);
      dimuon_mu1pt.push_back(vmu1.Pt()); dimuon_mu1eta.push_back(vmu1.Eta());
      dimuon_mu1phi.push_back(vmu1.Phi()); dimuon_mu2pt.push_back(vmu2.Pt()); 
      dimuon_mu2eta.push_back(vmu2.Eta()); dimuon_mu2phi.push_back(vmu2.Phi());
      dimuon_pt.push_back((vmu1+vmu2).Pt()); dimuon_eta.push_back((vmu1+vmu2).Eta());
      dimuon_phi.push_back((vmu1+vmu2).Phi()); dimuon_mass.push_back((vmu1+vmu2).M());
      GlobalPoint Dispbeamspot(-1*((beam_x-dimuvtx.position().x())+(dimuvtx.position().z()-beam_z)* theBeamSpot->dxdz()),-1*((beam_y-dimuvtx.position().y())+ (dimuvtx.position().z()-beam_z) * theBeamSpot->dydz()), 0);
      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
      math::XYZVector pperp(dimuvtx.position().x(),dimuvtx.position().y(),0);
      float cos(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
      dimuon_cos2D.push_back(cos);
      dimuon_lxy.push_back(Dispbeamspot.perp());
      dimuon_elxy.push_back(TMath::Sqrt(dimuvtx.positionError().rerr(Dispbeamspot)));
//      cout<<std::distance(dimucands.begin(),mu1)<<endl;
      dimuon_mu1id.push_back(dimuindex[std::distance(dimucands.begin(),mu1)]);
      dimuon_mu2id.push_back(dimuindex[std::distance(dimucands.begin(),mu2)]);
      ndimuons++;
   
   }
  }
 }
 if (MuMuVtx && SkipEvtNoMuVtx && !MuTrkVtx && dimuon_pt.size()==0) return;

 if (MuTrkVtx){
  for (std::vector<reco::Muon>::iterator mu1=dimucands.begin(); mu1!=dimucands.end(); ++mu1){
   for (std::vector<reco::Track>::iterator trk=trkcands.begin(); trk!=trkcands.end(); ++trk){ 
    if(mu1->charge()==trk->charge()) continue;
    if ( fabs(mu1->vz()-trk->vz())>DzMuTrk && DzMuTrk>0 ) continue;
    TLorentzVector vmu1,vmu2; 
    vmu1.SetPtEtaPhiM(mu1->pt(),mu1->eta(),mu1->phi(),0.105);
    vmu2.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(),0.105);
    if ( (vmu1+vmu2).M()<McutMin || (vmu1+vmu2).M()>McutMax) continue;
    std::vector<reco::TransientTrack> trkmutrk; trkmutrk.reserve(2);
    trkmutrk.emplace_back(reco::TransientTrack(*(mu1->bestTrack()),&(*bFieldHandle)));
    trkmutrk.emplace_back(reco::TransientTrack(*(trk),&(*bFieldHandle)));
    KalmanVertexFitter vtxFitter(true);
    TransientVertex trkmuvtx=vtxFitter.vertex(trkmutrk);
    if(!trkmuvtx.isValid()) continue; 
    if (ChiSquaredProbability(trkmuvtx.totalChiSquared(), trkmuvtx.degreesOfFreedom())<ProbCut) continue;
    GlobalPoint Dispbeamspot(-1*((beam_x-trkmuvtx.position().x())+(trkmuvtx.position().z()-beam_z)* theBeamSpot->dxdz()),-1*((beam_y-trkmuvtx.position().y())+ (trkmuvtx.position().z()-beam_z) * theBeamSpot->dydz()), 0);
    math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
    math::XYZVector pperp(trkmuvtx.position().x(),trkmuvtx.position().y(),0);
    float cos(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
    if (cos<CosCut) continue;
    trkmuon_mass_unfit.push_back((vmu1+vmu2).M());
    trkmuon_x.push_back(trkmuvtx.position().x());
    trkmuon_y.push_back(trkmuvtx.position().y());
    trkmuon_z.push_back(trkmuvtx.position().z());
    trkmuon_prob.push_back(ChiSquaredProbability(trkmuvtx.totalChiSquared(), trkmuvtx.degreesOfFreedom()));
    std::vector<reco::TransientTrack> refited=trkmuvtx.refittedTracks();
    vmu1.SetPtEtaPhiM(refited[0].track().pt(),refited[0].track().eta(),refited[0].track().phi(),0.105);
    vmu2.SetPtEtaPhiM(refited[1].track().pt(),refited[1].track().eta(),refited[1].track().phi(),0.105);
    trkmuon_mupt.push_back(vmu1.Pt()); trkmuon_mueta.push_back(vmu1.Eta());
    trkmuon_muphi.push_back(vmu1.Phi()); trkmuon_trkpt.push_back(vmu2.Pt()); 
    trkmuon_trketa.push_back(vmu2.Eta()); trkmuon_trkphi.push_back(vmu2.Phi());
    trkmuon_pt.push_back((vmu1+vmu2).Pt()); trkmuon_eta.push_back((vmu1+vmu2).Eta());
    trkmuon_phi.push_back((vmu1+vmu2).Phi()); trkmuon_mass.push_back((vmu1+vmu2).M());
    
    trkmuon_cos2D.push_back(cos);
    trkmuon_lxy.push_back(Dispbeamspot.perp());
    trkmuon_elxy.push_back(TMath::Sqrt(trkmuvtx.positionError().rerr(Dispbeamspot)));
    trkmuon_muid.push_back(dimuindex[std::distance(dimucands.begin(),mu1)]);
    trkmuon_trkid.push_back(trkindex[std::distance(trkcands.begin(),trk)]);
    ntrkmuons++;
  }
 }
}
if (MuTrkVtx && SkipEvtNoTrkVtx && !MuMuVtx && trkmuon_pt.size()==0 ) return;
if (MuTrkVtx && MuMuVtx && trkmuon_pt.size()==0 && dimuon_pt.size()==0 && SkipEvtNoTrkVtx && SkipEvtNoMuVtx) return;
  t1->Fill();


}


// ------------ method called once each job just before starting event loop  ------------
template<typename T1>
void 
MuonAODAnalyzer<T1>::beginJob()
{
 t1=fs->make<TTree>("mytree","mytree");
 t1->Branch("event",&event); t1->Branch("run_number",&run_number);
 t1->Branch("ls",&ls);
 t1->Branch("pvertex_x",&pvertex_x); t1->Branch("pvertex_y",&pvertex_y);
 t1->Branch("pvertex_z",&pvertex_z);
 //gen
 t1->Branch("genmuon_pt",&genmuon_pt); t1->Branch("genmuon_eta",&genmuon_eta);
 t1->Branch("genmuon_phi",&genmuon_phi); t1->Branch("genmuon_charge",&genmuon_charge);
 t1->Branch("genmuon_momId",&genmuon_momId); t1->Branch("genmuon_grandmomId",&genmuon_grandmomId);
 t1->Branch("genmuon_grandgrandmomId",&genmuon_grandgrandmomId);
 //trigger
 t1->Branch("hltpaths",&hltpaths);  
 t1->Branch("l1muon_pt",&l1muon_pt); t1->Branch("l1muon_eta",&l1muon_eta);
 t1->Branch("l1muon_phi",&l1muon_phi);
 t1->Branch("hltmuon_pt",&hltmuon_pt); t1->Branch("hltmuon_eta",&hltmuon_eta);
 t1->Branch("hltmuon_phi",&hltmuon_phi); t1->Branch("hltmupn_pathId",&hltmuon_pathId);
 t1->Branch("hlt_firedMuindex",&hlt_firedMuindex);
 t1->Branch("hlt_firedMuDr",&hlt_firedMuDr);
 //reco
 t1->Branch("nmuons",&nmuons); 
 t1->Branch("muon_pt",&muon_pt); t1->Branch("muon_eta",&muon_eta);
 t1->Branch("muon_phi",&muon_phi); t1->Branch("muon_charge",&muon_charge); 
 t1->Branch("muon_vx",&muon_vx); t1->Branch("muon_vy",&muon_vy);
 t1->Branch("muon_vz",&muon_vz); t1->Branch("muon_iso",&muon_iso);  
 t1->Branch("muon_soft",&muon_soft); t1->Branch("muon_loose",&muon_loose);
 t1->Branch("muon_medium",&muon_medium); t1->Branch("muon_tight",&muon_tight);
 t1->Branch("muon_dxy",&muon_dxy); t1->Branch("muon_edxy",&muon_edxy);
 t1->Branch("muon_dz",&muon_dz); t1->Branch("muon_edz",&muon_edz);
 t1->Branch("muon_global",&muon_global); t1->Branch("muon_pf",&muon_pf);
 t1->Branch("muon_normChi2",&muon_normChi2); t1->Branch("muon_nvalidhits",&muon_nvalidhits);
 t1->Branch("muon_nmatchstations",&muon_nmatchstations);
 t1->Branch("muon_npixelhits",&muon_npixelhits); t1->Branch("muon_ntrklayers",&muon_trklayers);
 t1->Branch("muon_vtxcand",&muon_vtxcand);
 //track
 t1->Branch("ntracks",&ntracks);
 t1->Branch("track_pt",&track_pt); t1->Branch("track_eta",&track_eta);
 t1->Branch("track_phi",&track_phi); t1->Branch("track_charge",&track_charge); 
 t1->Branch("track_vx",&track_vx); t1->Branch("track_vy",&track_vy);
 t1->Branch("track_vz",&track_vz); t1->Branch("track_dxy",&track_dxy);  
 t1->Branch("track_edxy",&track_edxy); t1->Branch("track_dz",&track_dz);  
 t1->Branch("track_edz",&track_edz); t1->Branch("track_trklayers",&track_trklayers);
 t1->Branch("track_pxllayers",&track_pxllayers); t1->Branch("track_normchi2",&track_normchi2);
 t1->Branch("track_vtxcand",&track_vtxcand);
 //dimuon
 t1->Branch("ndimuons",&ndimuons);
 t1->Branch("dimuon_pt",&dimuon_pt); t1->Branch("dimuon_eta",&dimuon_eta);
 t1->Branch("dimuon_phi",&dimuon_phi);  t1->Branch("dimuon_lxy",&dimuon_lxy);
 t1->Branch("dimuon_elxy",&dimuon_elxy); t1->Branch("dimuon_cos2D",&dimuon_cos2D);
 t1->Branch("dimuon_x",&dimuon_x); t1->Branch("dimuon_y",&dimuon_y);
 t1->Branch("dimuon_z",&dimuon_z); t1->Branch("dimuon_prob",&dimuon_prob);
 t1->Branch("dimuon_mass",&dimuon_mass); t1->Branch("dimuon_mass_unfit",&dimuon_mass_unfit);
 t1->Branch("dimuon_mu1pt",&dimuon_mu1pt); t1->Branch("dimuon_mu1eta",&dimuon_mu1eta);
 t1->Branch("dimuon_mu1phi",&dimuon_mu1phi); t1->Branch("dimuon_mu2pt",&dimuon_mu2pt);
 t1->Branch("dimuon_mu2eta",&dimuon_mu2eta); t1->Branch("dimuon_mu2phi",&dimuon_mu2phi);
 t1->Branch("dimuon_mu1id",&dimuon_mu1id); t1->Branch("dimuon_mu2id",&dimuon_mu2id);
 //track muon
 t1->Branch("ntrkmuons",&ntrkmuons);
 t1->Branch("trkmuon_pt",&trkmuon_pt); t1->Branch("trkmuon_eta",&trkmuon_eta);
 t1->Branch("trkmuon_phi",&trkmuon_phi);  t1->Branch("trkmuon_lxy",&trkmuon_lxy);
 t1->Branch("trkmuon_elxy",&trkmuon_elxy); t1->Branch("trkmuon_cos2D",&trkmuon_cos2D);
 t1->Branch("trkmuon_x",&trkmuon_x); t1->Branch("trkmuon_y",&trkmuon_y);
 t1->Branch("trkmuon_z",&trkmuon_z); t1->Branch("trkmuon_prob",&trkmuon_prob);
 t1->Branch("trkmuon_mass",&trkmuon_mass); t1->Branch("trkmuon_mass_unfit",&trkmuon_mass_unfit);
 t1->Branch("trkmuon_mupt",&trkmuon_mupt); t1->Branch("trkmuon_mueta",&trkmuon_mueta);
 t1->Branch("trkmuon_muphi",&trkmuon_muphi); t1->Branch("trkmuon_trkpt",&trkmuon_trkpt);
 t1->Branch("trkmuon_trketa",&trkmuon_trketa); t1->Branch("trkmuon_trkphi",&trkmuon_trkphi);
 t1->Branch("trkmuon_muid",&trkmuon_muid); t1->Branch("trkmuon_trkid",&trkmuon_trkid);
}

// ------------ method called once each job just after ending the event loop  ------------
template<typename T1>
void 
MuonAODAnalyzer<T1>::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<typename T1>
void
MuonAODAnalyzer<T1>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


///////////////////////
  
//define this as a plug-in
typedef MuonAODAnalyzer<reco::RecoEcalCandidate> MuonAODAnalyzerb;
DEFINE_FWK_MODULE(MuonAODAnalyzerb);


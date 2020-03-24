import FWCore.ParameterSet.Config as cms
#quick config


IsData=False
Nentries=5000;  output="output_check.root"; nlog=1000; #Debug='MIN:MIN:MIN-MAX:MAX:MAX'
Dataset='SingleMuon' #options: SingleMuon or BParking or 5TeV
DrMuTrgMatch=dict(SkipIfNoMatch=False,Cone=0.05); 
SkipEvtWithNoFire=False; 
SaveTrk=dict(PassSelection=False,MatchedToGen=True)

#Objects
#determine muon/track cands for tag&probe. The muons col. in ntuple will not "respect" those cuts (ie all muons saved). Cuts *just* limiting muons in t&p module
Preselection=dict(Gen=False,Trg=False,Cut=False); #select only one
# 0 use directly gen mu and take close muons tracks
GenPresel=dict(MomPDG=531,GenRecoCone=0.1)
# 1 use directly trg muon - no extra options needed, obviously only in 1 muon
# 2 Apply preselection cuts
MuPresel=dict(Pt=20.0,DzTrg=1.0,DrTrg=0.4,Soft=False)
#cuts for track selection in case of trk-mu vtx - applied also in track_* in ntuple. If GenSelection is True on *top* of selection cuts only the tracks close to muon is provided as fitting candidates
TrkPresel=dict(Pt=1.0,DzTrg=1.0) #gen selection is excluded from those cuts


MinM=2.7; MaxM=3.3;

MuMu=dict(PtMin=0,DzMuMu=10,Vtx=False); 
MuTrk=dict(PtMin=-3,Prob=-0.0001,Cos=-10.0,DzMuTrk=10.0,Vtx=False);
SkipEvtWithNoMuVtx=False; SkipEvtWithNoTrkVtx=False;
File=[
#'/store/data/Run2018C/SingleMuon/AOD/PromptReco-v3/000/319/993/00000/0603799C-778E-E811-8792-FA163EAD4CB1.root'
#'/store/mc/RunIIAutumn18DR/BuToKJpsi_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/810000/F647981D-B4EE-E74E-9CC6-A7224B0FB8DD.root'
#'/store/mc/RunIIAutumn18DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/80002/FF9AF238-78B6-CF48-BC7C-05025D85A45C.root'
#'/store/himc/RunIIpp5Spring18DR/JPsiMM_TuneCUETP8M1_5p02TeV_pythia8/AODSIM/94X_mc2017_realistic_forppRef5TeV-v2/60000/F6DB806A-CF46-E911-A738-B083FED406AC.root'
#'/store/mc/RunIIAutumn18DRPremix/BsToMuMuPhi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/102X_upgrade2018_realistic_v15-v1/120000/ADA9C6F9-EF32-5747-B816-8F850FA1922B.root'
#'/store/data/Run2017G/SingleMuon/AOD/17Nov2017-v1/90002/183672E6-A42E-E811-92C1-0025904AC2CC.root'
#'/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/60004/FB123080-071C-F64D-BAFD-F2F292F7FC64.root'
  "/store/user/eliza/MCjpsiMuMuTNP2018.root"
]

############
Path=[]; Filter=[];   

if Dataset=="SingleMuon":
  #Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu19_v","HLT_Mu20_v","HLT_IsoMu20_v","HLT_IsoMu24_v","HLT_Mu50"]; 
  #Filter=["hltL3fL1sMu5L1f0L2f5L3Filtered8","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered19","hltL3fL1sMu18L1f0L2f10QL3Filtered20Q","hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p07","hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"]
  Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu20_v"];
  Filter=["hltL3fL1sMu5L1f0L2f5L3Filtered8","hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17","hltL3fL1sMu18L1f0L2f10QL3Filtered20Q"]

if Dataset=="BParking": 
  Path=["HLT_Mu9_IP6_part","HLT_Mu8p5_IP3p5","HLT_Mu10p5_IP3p5","HLT_Mu8_IP3"]
  Filter=["hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8p5Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered10p5Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q"]; 

if Dataset=="5TeV": 
  Path=["HLT_HIL2Mu7_v1","HLT_HIL2Mu12_v1","HLT_HIL2Mu15_v1","HLT_HIL2Mu20_v1", "HLT_HIL3Mu3_v1","HLT_HIL3Mu5_v1","HLT_HIL3Mu7_v1","HLT_HIL3Mu12_v1", "HLT_HIL3Mu15_v1","HLT_HIL3Mu20_v1","HLT_HIL2Mu3_NHitQ10_v1","HLT_HIL3Mu3_NHitQ10_v1","HLT_HIL2Mu5_NHitQ10_v1","HLT_HIL3Mu5_NHitQ10_v1"]
  Filter=["hltL2fL1sSingleMu3OR5L1f0L2Filtered7","hltL2fL1sSingleMu7L1f0L2Filtered12","hltL2fL1sSingleMu7L1f0L2Filtered15","hltL2fL1sSingleMu7L1f0L2Filtered20",
"hltL3fL1sSingleMu3L1f0L2f0L3Filtered3","hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered5","hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered7","hltL3fL1sSingleMu7L1f0L2f0L3Filtered12","hltL3fL1sSingleMu7L1f0L2f0L3Filtered15","hltL3fL1sSingleMu7L1f0L2f0L3Filtered20","hltL2fL1sSingleMu3L1f0L2NHitQ10L2Filtered3","hltL3fL1sSingleMu3L1f0L2f0L3NHitQ10L3Filtered3","hltL2fL1sSingleMu3OR5L1f0L2NHitQ10L2Filtered5","hltL3fL1sSingleMu3OR5L1f0L2f0L3NHitQ10L3Filtered5"]; 


globaltag='102X_upgrade2018_realistic_v15' 
if IsData:
   print "We have established we Run on data using HLT",Path
   globaltag='101X_dataRun2_Prompt_v11'     
else:
   print "We have established we Run on MC"

print "Run parameters ",globaltag,""

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = nlog

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(Nentries) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    File ),
   secondaryFileNames=cms.untracked.vstring(),
#   eventsToProcess=cms.untracked.VEventRange("316187:827:MIN-316187:827:MAX"),
   inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_ctppsPixelClusters_*_*')

)

process.demo = cms.EDAnalyzer('MuonAODAnalyzerb',
                              beamSpot=cms.InputTag('offlineBeamSpot'),
                              vertices=cms.InputTag("offlinePrimaryVertices"),
                              triggerresults=cms.InputTag("TriggerResults::HLT"),
                              triggerobjects=cms.InputTag('hltTriggerSummaryAOD'), 
                              muons=cms.InputTag("muons"),
                              l1muons=cms.InputTag("gmtStage2Digis","Muon","RECO"),
                              tracks=cms.InputTag("generalTracks"),
         
                              gen = cms.InputTag("genParticles"),
                              HLTPaths=cms.vstring(Path),
                              HLTFilters=cms.vstring(Filter),
                               RunParameters = cms.PSet(
      IsData= cms.bool(IsData),
      SkipEvtWithNoFire=cms.bool(SkipEvtWithNoFire),
      SaveTracks=cms.bool(SaveTrk["PassSelection"]),
      TrgMuMatchCone=cms.double(DrMuTrgMatch["Cone"]),
      SkipEvtNoHLTmuMatch=cms.bool(DrMuTrgMatch["SkipIfNoMatch"]),
      CutSelection=cms.bool(Preselection['Cut']),
      TrgSelection=cms.bool(Preselection['Trg']),
      GenSelection=cms.bool(Preselection['Gen']),
      #mu cuts only for cands
      MuPtCut=cms.double(MuPresel['Pt']),
      DzMuTrgCut=cms.double(MuPresel['DzTrg']),
      DrMuTrgCut=cms.double(MuPresel['DrTrg']),
      SoftQCut=cms.bool(MuPresel['Soft']),
      #track preselection 
      TrkPtCut=cms.double(TrkPresel['Pt']),
      DzTrkTrg=cms.double(TrkPresel['DzTrg']),
      #gen preselection
      MomId=cms.double(GenPresel['MomPDG']),
      GenRecoMatchCone=cms.double(GenPresel['GenRecoCone']),
      SaveMatchGenTrkOnly=cms.bool(SaveTrk["MatchedToGen"]),
      # mumu vertex
      MuMuVtx=cms.bool(MuMu["Vtx"]),
      DzMuMu=cms.double(MuMu['DzMuMu']),
      minDiMuPtCut=cms.double(MuMu["PtMin"]),
      # mutrk vertex
      MuTrkVtx=cms.bool(MuTrk["Vtx"]),
      DzMuTrk=cms.double(MuTrk['DzMuTrk']),
      minMuTrkPtCut=cms.double(MuTrk["PtMin"]),
      MuTrkProbCut=cms.double(MuTrk["Prob"]),
      MuTrkCosCut=cms.double(MuTrk["Cos"]),
      # common 
      MresonanceMin=cms.double(MinM),
      MresonanceMax=cms.double(MaxM),
      SkipEvtNoMuVtx=cms.bool(SkipEvtWithNoMuVtx),
      SkipEvtNoTrkVtx=cms.bool(SkipEvtWithNoTrkVtx),
     
  ),
)



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(output)
                                   )
process.fevt = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(#"drop *",
    ),
    fileName = cms.untracked.string("edm_output.root"))

process.p = cms.Path( process.demo )
   

import FWCore.ParameterSet.Config as cms
#quick config


IsData=True
Nentries=-1;  output="output_mc.root"; nlog=100; #Debug='MIN:MIN:MIN-MAX:MAX:MAX'
#Single Muon paths
#Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu19_v","HLT_Mu20_v","HLT_IsoMu20_v"]; 
#Bphysics 
Path=["HLT_Mu9_IP6","HLT_Mu7_IP4"]
#Path=["HLT_Mu50","HLT_TrkMu50"]; 
DrTrgCone=dict(SkipIfNoMatch=True,Max=0.05); 
SaveTrk=True
MinM=2.6; MaxM=3.3; MuPtMin=8.0; TrkPtMin=1.5; DzMuTrk=10; DzMuMu=1000; DzTrkTrg=10;
UseOnlyTrgForTrkMu=True
MuMu=dict(PtMin=0,Vtx=False); MuTrk=dict(PtMin=5,Vtx=True);
SkipEvtWithNoFire=True; SkipEvtWithNoMuVtx=False; SkipEvtWithNoTrkVtx=True;
HighPtMu=dict(Thr=500,SkipLowPtEvt=False);
File=[
#'/store/mc/RunIIAutumn18MiniAOD/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/120000/078DB2B1-40DD-634D-A3CF-D2E377CAFA48.root'
#'/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/80000/FA6CC04A-6C3C-6C4D-89BC-E2B924290010.root'
#'/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/80000/FA6CC04A-6C3C-6C4D-89BC-E2B924290010.root'
#'/store/data/Run2018D/ZeroBias9/MINIAOD/PromptReco-v2/000/324/578/00000/FA330CBB-3E49-2340-98F5-C46002C02B13.root'
#'/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/00000/0109D8A3-FF62-6245-9745-675EE9FD6243.root'
#'/store/data/Run2018B/SingleMuon/MINIAOD/17Sep2018-v1/60000/FEC8A983-6AFB-2A44-AEE9-C2DC88E7C8F1.root'
#'/store/data/Run2018A/ParkingBPH6/MINIAOD/14May2018-v1/10000/8CBC6070-D558-E811-8D6F-A4BF0112BC28.root'
'/store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/F48C8243-CA59-E811-B3F9-FA163EBB80A0.root'
#data
#'/store/data/Run2017E/SingleMuon/MINIAOD/17Nov2017-v1/60000/F2A4FC7F-57DD-E711-BBB9-02163E019E29.root'
#'/store/data/Run2017C/SingleMuon/MINIAOD/17Nov2017-v1/40000/04BC4AEA-6BD8-E711-B000-02163E01A450.root'
#'/store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/110000/BBEEDE43-EE76-9F4D-A9C3-CFD1D2E58753.root'
#'/store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/110000/47699061-1903-F847-82ED-5E164D41DF96.root'
]

############
   

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
   eventsToProcess=cms.untracked.VEventRange("316187:827:MIN-316187:827:MAX"),
   inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_ctppsPixelClusters_*_*')

)

process.demo = cms.EDAnalyzer('MuonAnalyzerb',
                              beamSpot=cms.InputTag('offlineBeamSpot'),
                              vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
                              triggerresults=cms.InputTag("TriggerResults::HLT"),
                              triggerobjects=cms.InputTag('slimmedPatTrigger'), 
                              muons=cms.InputTag("slimmedMuons"),
                              l1muons=cms.InputTag("gmtStage2Digis","Muon","RECO"),
                              PFCands=cms.InputTag("packedPFCandidates"),
                              losttracks=cms.InputTag("lostTracks"),
                              genmuons = cms.InputTag("packedGenParticles"),
                              HLTPaths=cms.vstring(Path),

                               RunParameters = cms.PSet(
      Data= cms.bool(IsData),SkipEvtWithNoFire=cms.bool(SkipEvtWithNoFire),
      SaveTracks=cms.bool(SaveTrk),
      maxDrMuMatch=cms.double(DrTrgCone["Max"]),
      SkipEvtNoHLTmuMatch=cms.bool(DrTrgCone["SkipIfNoMatch"]),
      minMuPtCut=cms.double(MuPtMin),minTrkPtCut=cms.double(TrkPtMin),
      DzMuTrk=cms.double(DzMuTrk),DzMuMu=cms.double(DzMuMu),DzTrkTrg=cms.double(DzTrkTrg),
      MuMuVtx=cms.bool(MuMu["Vtx"]),MuTrkVtx=cms.bool(MuTrk["Vtx"]),
      McutMin=cms.double(MinM),McutMax=cms.double(MaxM),
      minDiMuPtCut=cms.double(MuMu["PtMin"]),
      minMuTrkPtCut=cms.double(MuTrk["PtMin"]),
      SkipEvtNoMuVtx=cms.bool(SkipEvtWithNoMuVtx),
      SkipEvtNoTrkVtx=cms.bool(SkipEvtWithNoTrkVtx),
      UseOnlyTrgForTrkMu=cms.bool(UseOnlyTrgForTrkMu),
      HighPtThr=cms.double(HighPtMu["Thr"]),
      SkipEvtNoHighPt=cms.bool(HighPtMu["SkipLowPtEvt"])
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
   

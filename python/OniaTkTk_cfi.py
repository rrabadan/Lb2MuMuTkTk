import FWCore.ParameterSet.Config as cms


PsiTTProducer = cms.EDProducer('OniaTrakTrakProducer',
    Onia = cms.InputTag('onia2MuMuPAT'),
    Trak = cms.InputTag('selectedPatTrackCands'),#clean
    dEdx1Tag = cms.InputTag('dedxHarmonic2'),
    dEdx2Tag = cms.InputTag('dedxDiscrimASmi'),
    OniaMassCuts = cms.vdouble(2.946916,3.246916),  
    TrakTrakMassCuts = cms.vdouble(2.,999),
    OniaTrakTrakMassCuts = cms.vdouble(2.,999),            # b-hadron mass window
    MassTraks    = cms.vdouble(0.938272,0.493677),         # traks masses
    OnlyBest  = cms.bool(False)    
)

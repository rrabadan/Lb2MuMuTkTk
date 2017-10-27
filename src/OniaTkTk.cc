#include "Lb2MuMuTrakTrak/interface/OniaTrakTrakProducer.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Ref.h"

OniaTrakTrakProducer::OniaTrakTrakProducer(const edm::ParameterSet& ps):
OniaCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("Onia"))),
TrakCollection_(consumes<std::vector<pat::GenericParticle>>(ps.getParameter<edm::InputTag>("Trak"))),
m_dEdx1Tag_(consumes<edm::ValueMap<reco::DeDxData> >(ps.getParameter< edm::InputTag >("dEdx1Tag")));
m_dEdx2Tag_(consumes<edm::ValueMap<reco::DeDxData> >(ps.getParameter< edm::InputTag >("dEdx2Tag")));
OniaMassCuts_(ps.getParameter<std::vector<double>>("OniaMassCuts")),
TrakTrakMassCuts_(ps.getParameter<std::vector<double>>("TrakTrakMassCuts")),
OniaTrakTrakMassCuts_(ps.getParameter<std::vector<double>>("OniaTrakTrakMassCuts")),
MassTraks_(ps.getParameter<std::vector<double>>("MassTraks")),
OnlyBest_(ps.getParameter<bool>("OnlyBest"))
{
    produces<pat::CompositeCandidateCollection>("Candidates");
    candidates = 0;
    nevents = 0;
    nonia = 0;
    nreco = 0;
}

void OniaTrakTrakProducer::produce(edm::Event& event, const edm::EventSetup& esetup){
    
    std::auto_ptr<pat::CompositeCandidateCollection> OniaTTCandColl(new pat::CompositeCandidateCollection);
    
    edm::Handle<pat::CompositeCandidateCollection> onia;
    event.getByToken(OniaCollection_,onia);
    
    edm::Handle<std::vector<pat::GenericParticle> > trak;
    event.getByToken(TrakCollection_,trak);
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdx1Handle;
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdx2Handle;
    event.getByToken(dEdx1Tag_,dEdx1Handle);
    event.getByToken(dEdx2Tag_,dEdx2Handle);
    std::vector<const edm::ValueMap<reco::DeDxData> *> v_dEdx;
    v_dEdx.push_back(dEdx1Handle.product());
    v_dEdx.push_back(dEdx2Handle.product());
    
    uint ncombo = 0;
    float OniaMassMax_ = OniaMassCuts_[1];
    float OniaMassMin_ = OniaMassCuts_[0];
    float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
    float TrakTrakMassMin_ = TrakTrakMassCuts_[0];
    float OniaTrakTrakMassMax_ = OniaTrakTrakMassCuts_[1];
    float OniaTrakTrakMassMin_ = OniaTrakTrakMassCuts_[0];
    
    // Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
    for (pat::CompositeCandidateCollection::const_iterator oniaCand = onia->begin(); oniaCand != onia->end(); ++oniaCand){
        if ( oniaCand->mass() < OniaMassMax_  && oniaCand->mass() > OniaMassMin_ ) {
            const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon1"));
            const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon2"));
            if ( pmu1->pt() < 2.9) continue;
            if ( pmu2->pt() < 2.9) continue;
            for (std::vector<pat::GenericParticle>::const_iterator trakCand = trak->begin(), trakend=trak->end(); trakCand!= trakend; ++trakCand){
                if ( IsTheSame(*trakCand,*pmu1) || IsTheSame(*trakCand,*pmu2) ) continue;
                if ( trakCand->pt() < 0.9) continue;
                if ( trakCand->track()->chi2()>30) continue;
                // loop over second track candidate, negative charge
                for (std::vector<pat::GenericParticle>::const_iterator trakCand2 = trakCand+1; trakCand2!= trakend; ++trakCand2){
                    if (trakCand2 == trakCand) continue;
                    if ( trakCand2->charge() == trakCand->charge()) continue;
                    //if ( IsTheSame(*trakCand2,*pmu1) || IsTheSame(*trakCand2,*pmu2) || trakCand2->charge() > 0 ) continue;
                    if ( IsTheSame(*trakCand2,*pmu1) || IsTheSame(*trakCand2,*pmu2) ) continue;
                    if ( trakCand2->pt() < 0.9) continue;
                    if ( trakCand2->track()->chi2()>30) continue;
                    pat::CompositeCandidate TTCand = makeTTCandidate(*trakCand, *trakCand2, MassTraks_);
                    if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {
                        pat::CompositeCandidate OniaTTCand = makeOniaTTCandidate(*oniaCand, *&TTCand);
                
                        if ( OniaTTCand.mass() < OniaTrakTrakMassMax_ && OniaTTCand.mass() > OniaTrakTrakMassMin_) {
                            //OniaTTCandColl->push_back(OniaTTCand);
                            pat::CompositeCandidate TTCand2 = makeTTCandidate(*trakCand2, *trakCand, MassTraks_);
                            pat::CompositeCandidate OniaTTCand2 = makeOniaTTCandidate(*oniaCand, *&TTCand2);
                            
                            
                            float x2_1 = trakCand->track()->normalizedChi2();
                            float x2_2 = trakCand2->track()->normalizedChi2();
                            
                            float dEdx_1[2] = {0.0,0.0};
                            float dEdx_2[2] = {0.0,0.0};
                            for (unsigned int i=0; i<v_dEdx.size(); i++) {
                                const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[i]);
                                const reco::DeDxData& dedx_1 = dEdxTrack[trakCand->track()];
                                const reco::DeDxData& dedx_2 = dEdxTrack[trakCand2->track()];
                                dEdx_1[i]= dedx_1.dEdx();
                                dEdx_2[i]= dedx_2.dEdx();
                            }
                            
                            OniaTTCand.addUserFloat("p_chi2", x2_1);
                            OniaTTCand.addUserFloat("p_dEdxH", dEdx_1[0]);
                            OniaTTCand.addUserFloat("p_dEdxD", dEdx_1[1]);
                            
                            
                            OniaTTCand2.addUserFloat("p_chi2", x2_2);
                            OniaTTCand2.addUserFloat("k_dEdxH", dEdx_1[0]);
                            OniaTTCand2.addUserFloat("k_dEdxD", dEdx_1[1]);
                            //std::cout << "onia: " << oniaCand->mass() << " cand1: " << OniaTTCand.mass() << " cand2: " << OniaTTCand2.mass() << "nomChi2: " << x2_1 <<std::endl;
                            
                            OniaTTCandColl->push_back(OniaTTCand);
                            OniaTTCandColl->push_back(OniaTTCand2);
                        }
                    }
                } // loop over second track
            }   // loop on track candidates
        }
        if (OnlyBest_) break;
    }
    event.put(OniaTTCandColl,"Candidates");
    nevents++;
}

void OniaTrakTrakProducer::endJob(){
    
}

bool OniaTrakTrakProducer::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
    double DeltaEta = fabs(mu.eta()-tk.eta());
    double DeltaP   = fabs(mu.p()-tk.p());
    if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
    return false;
}

const pat::CompositeCandidate OniaTrakTrakProducer::makeOniaTTCandidate(
                                                                        const pat::CompositeCandidate& onia,
                                                                        const pat::CompositeCandidate& tt
                                                                        ){
    
    pat::CompositeCandidate OniaTCand;
    OniaTCand.addDaughter(onia,"onia");
    OniaTCand.addDaughter(tt,"ditrak");
    OniaTCand.setVertex(onia.vertex());
    OniaTCand.setCharge(tt.charge());
    
    reco::Candidate::LorentzVector vOniaT = onia.p4() + tt.p4();
    OniaTCand.setP4(vOniaT);
    
    return OniaTCand;
    
}

const pat::CompositeCandidate OniaTrakTrakProducer::makeTTCandidate(
                                                                    const pat::GenericParticle& trak1,
                                                                    const pat::GenericParticle& trak2, const std::vector<double>& masstraks_
                                                                    ){
    
    pat::CompositeCandidate TTCand;
    TTCand.addDaughter(trak1,"trak1");
    TTCand.addDaughter(trak2,"trak2");
    TTCand.setCharge(trak1.charge()+trak2.charge());
    
    double m_trak1 = masstraks_[0];
    //std::cout << "---------------------------- mass0 " <<  masstraks_[0] << std::endl;
    math::XYZVector mom_trak1 = trak1.momentum();
    double e_trak1 = sqrt(m_trak1*m_trak1 + mom_trak1.Mag2());
    math::XYZTLorentzVector p4_trak1 = math::XYZTLorentzVector(mom_trak1.X(),mom_trak1.Y(),mom_trak1.Z(),e_trak1);
    double m_trak2 = masstraks_[1];
    math::XYZVector mom_trak2 = trak2.momentum();
    double e_trak2 = sqrt(m_trak2*m_trak2 + mom_trak2.Mag2());
    math::XYZTLorentzVector p4_trak2 = math::XYZTLorentzVector(mom_trak2.X(),mom_trak2.Y(),mom_trak2.Z(),e_trak2);
    reco::Candidate::LorentzVector vTT = p4_trak1 + p4_trak2;
    TTCand.setP4(vTT);
    
    return TTCand;
}

reco::Candidate::LorentzVector OniaTrakTrakProducer::convertVector(const math::XYZTLorentzVectorF& v){
    
    return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(OniaTrakTrakProducer);

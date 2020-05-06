//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// helper functions


#ifndef MuonAnalysis_MuonAnalyzer_plugins_helper
#define MuonAnalysis_MuonAnalyzer_plugins_helper

const float MU_MASS=0.10565837;
inline float DimuonMass(float mu1pt,float mu1eta,float mu1phi,
                        float mu2pt,float mu2eta,float mu2phi)
       {
         math::PtEtaPhiMLorentzVector mu1(mu1pt, mu1eta, mu1phi, MU_MASS);
         math::PtEtaPhiMLorentzVector mu2(mu2pt, mu2eta, mu2phi, MU_MASS);
         return (mu1+mu2).mass();              
       }


#endif

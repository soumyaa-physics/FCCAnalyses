#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/VertexFitterSimple.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/ReconstructedTrack.h"

namespace FCCAnalyses{

namespace ReconstructedParticle2Track{

/// Kink Finder
ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> findKink_candidate(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced,
    const ROOT::VecOps::RVec<edm4hep::TrackState> & fullTracks 
){
   // full tracks is needed to get the index
   // return the vertex data of the kink candidate, if any. If multiple candidates, return all of them. If none, return empty vector.
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> vertexCandidates;
    
    auto primary_indices = ReconstructedTrack::get_indices(primary, fullTracks);
    auto displaced_indices = ReconstructedTrack::get_indices(displaced, fullTracks);
    // useful for background
    if (displaced.empty()) return vertexCandidates;

    auto primMom_mag = getRP2TRK_mom(reco, primary); // this gives the magnitude of the momentum
    auto dispMom_mag = getRP2TRK_mom(reco, displaced);

    // moved this out of the loop, as it only needs to be done once
    std::vector<size_t> dispIdx(displaced.size()); //size of displaced tracks
    std::iota(dispIdx.begin(), dispIdx.end(), 0); // initialize with 0,1,2,
    // sort indices of displaced tracks by momentum
    std::sort(dispIdx.begin(), dispIdx.end(), [&](size_t a, size_t b) { return dispMom_mag[a] > dispMom_mag[b]; });
    size_t nDispToUse = std::min<size_t>(3, dispIdx.size()); // get the first three indices

    // loop over primary tracks
    for (size_t i = 0; i < primary.size(); ++i) {
        int trackIndex = primary_indices[i]; // get the correct index
        if (trackIndex < 0) continue; // if no valid index, skip
        const auto& tIn = primary[i];

        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) {
            for (int it = rp.tracks_begin; it < rp.tracks_end; ++it) {
                if (it == trackIndex) { // index of tracks and reco particles should match
                    rpIn = &rp;
                    break;
                }
            }
            if (rpIn) break; // if we found the associated reco particle, no need to continue looping
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        // auto primMom_current = primMom[i];
        bool kinkFound = false; 
        // doing the same for displaced tracks
        for (size_t k = 0; k < nDispToUse; ++k) {
            if (kinkFound) break; 
            size_t j = dispIdx[k];
            int dispTrackIndex = displaced_indices[j];
            if (dispTrackIndex < 0) continue; // if no valid index, skip

            const auto& tOut = displaced[j]; // current displaced track

            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                for (int it = rp.tracks_begin; it < rp.tracks_end; ++it) {
                    if (it == dispTrackIndex) {
                        rpOut = &rp;
                        break;
                    }
                }
                if (rpOut) break;
            }
            if (!rpOut) continue;
            // auto dispMom_current = dispMom[j]; // use the current displaced momentum

            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) < 0.05) continue;
            // vertex reconstruction
            ROOT::VecOps::RVec<edm4hep::TrackState> tracksToFit = { tIn, tOut };
            auto vtxObj = VertexFitterSimple::VertexFitter_Tk(2, tracksToFit); // flag 2 for SVs
            auto vtxData = VertexingUtils::get_VertexData(vtxObj);
            ROOT::VecOps::RVec<float> vertex = {
                vtxData.position.x,
                vtxData.position.y,
                vtxData.position.z
            };

            // removing duplicates- ig not needed when i only 2 leading displaced tracks, but just in case
            bool isDuplicate = false;

            for (auto& existing : vertexCandidates) {
                float dx = existing[0] - vertex[0];
                float dy = existing[1] - vertex[1];
                float dz = existing[2] - vertex[2];
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist < 0.1f) {
                    isDuplicate = true;
                    break;
                }
            }

            if (!isDuplicate)
            vertexCandidates.push_back(vertex);   
            kinkFound = true;           
        }
    }
    return vertexCandidates;
  }
/// Kink Finder- return vertex object
 ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> KinkCandidate_VertexObject(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced,
    const ROOT::VecOps::RVec<edm4hep::TrackState> & fullTracks
){
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertexCandidates;
    auto primary_indices = ReconstructedTrack::get_indices(primary, fullTracks);
    auto displaced_indices = ReconstructedTrack::get_indices(displaced, fullTracks);
    // useful for background
    if (displaced.empty()) return vertexCandidates;

    auto primMom_mag = getRP2TRK_mom(reco, primary); // this gives the magnitude of the momentum
    auto dispMom_mag = getRP2TRK_mom(reco, displaced);

    // moved this out of the loop, as it only needs to be done once
    std::vector<size_t> dispIdx(displaced.size()); //size of displaced tracks
    std::iota(dispIdx.begin(), dispIdx.end(), 0); // initialize with 0,1,2,
    // sort indices of displaced tracks by momentum
    std::sort(dispIdx.begin(), dispIdx.end(), [&](size_t a, size_t b) { return dispMom_mag[a] > dispMom_mag[b]; });
    size_t nDispToUse = std::min<size_t>(3, dispIdx.size()); // get the first three indices

    // loop over primary tracks
    for (size_t i = 0; i < primary.size(); ++i) {
        int trackIndex = primary_indices[i]; // get the correct index
        if (trackIndex < 0) continue; // if no valid index, skip
        const auto& tIn = primary[i];

        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) {
            for (int it = rp.tracks_begin; it < rp.tracks_end; ++it) {
                if (it == trackIndex) { // index of tracks and reco particles should match
                    rpIn = &rp;
                    break;
                }
            }
            if (rpIn) break; // if we found the associated reco particle, no need to continue looping
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        // auto primMom_current = primMom[i];
        bool kinkFound = false;  // flag per primary track

        // doing the same for displaced tracks

        for (size_t k = 0; k < nDispToUse; ++k) {
            if (kinkFound) break; 
            size_t j = dispIdx[k];
            int dispTrackIndex = displaced_indices[j];
            if (dispTrackIndex < 0) continue; // if no valid index, skip

            const auto& tOut = displaced[j]; // current displaced track

            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                for (int it = rp.tracks_begin; it < rp.tracks_end; ++it) {
                    if (it == dispTrackIndex) {
                        rpOut = &rp;
                        break;
                    }
                }
                if (rpOut) break;
            }
            if (!rpOut) continue;
            // auto dispMom_current = dispMom[j]; // use the current displaced momentum

            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) < 0.05) continue;
            // vertex reconstruction
            ROOT::VecOps::RVec<edm4hep::TrackState> tracksToFit = { tIn, tOut };
            auto vtxObj = VertexFitterSimple::VertexFitter_Tk(2, tracksToFit); // flag 2 for SVs
            auto vtxData = VertexingUtils::get_VertexData(vtxObj);
            ROOT::VecOps::RVec<float> vertex = {
                vtxData.position.x,
                vtxData.position.y,
                vtxData.position.z
            };
            vertexCandidates.push_back(vtxObj); 
            kinkFound = true;           
        }
    }
    return vertexCandidates;
  }

  ROOT::VecOps::RVec<float> 
  getRP2TRK_mom(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
                ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<float> result;
    for (auto & p: in) {
      if (p.tracks_begin<tracks.size())
        result.push_back(VertexingUtils::get_trackMom(tracks.at(p.tracks_begin)));
      else result.push_back(std::nan(""));
    }
    return result;
  }

  ROOT::VecOps::RVec<float> 
  getRP2TRK_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
                   ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<float> result;
    for (auto & p: in) {
      if (p.tracks_begin<tracks.size())
        result.push_back(p.charge);
      else result.push_back(std::nan(""));
    }
    return result;
  }

  ROOT::VecOps::RVec<float> getRP2TRK_Bz(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps, const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks) {
    const double c_light = 2.99792458e8;
    const double a = c_light * 1e3 * 1e-15; //[omega] = 1/mm
    ROOT::VecOps::RVec<float> out;

    for(auto & p: rps) {
      if(p.tracks_begin < tracks.size()) {
	double pt= sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y);
	double Bz= tracks.at(p.tracks_begin).omega / a * pt * std::copysign(1.0, p.charge);
	out.push_back(Bz);
      } else {
	out.push_back(-9.);
      }
    }
    return out;
  }

  float Bz(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps, const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks) {
    const double c_light =  2.99792458e8;// speed of light m/sec;
    const double a = c_light * 1e3 * 1e-15; //[omega] = 1/mm

    double Bz = -9;

    for(auto & p: rps) {
      if(p.tracks_begin < tracks.size()) {
        double pt= sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y);
        Bz= tracks.at(p.tracks_begin).omega / a * pt * std::copysign(1.0, p.charge);
      }
    }
    return Bz;
  }

  ROOT::VecOps::RVec<float> XPtoPar_dxy(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
					const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
					const TLorentzVector& V, // primary vertex
					const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e-9;

    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        float D0_wrt0 = tracks.at(rp.tracks_begin).D0;
        float Z0_wrt0 = tracks.at(rp.tracks_begin).Z0;
        float phi0_wrt0 = tracks.at(rp.tracks_begin).phi;

        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - V.Vect();
        //std::cout<<"vertex: "<<V.Vect().X()<<", "<<V.Vect().Y()<<", "<<V.Vect().Z()<<", "<<std::endl;
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);

        double a = - rp.charge * Bz * cSpeed;
        double pt = p.Pt();
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double D=-9;
        if (pt * pt - 2 * a * cross + a * a * r2 > 0) {
          double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
      	  if (pt < 10.0) D = (T - pt) / a;
                else D = (-2 * cross + a * r2) / (T + pt);
        }
        //std::cout<<"displ: "<<D<<std::endl;
	      out.push_back(D);

      } else {
	out.push_back(-9.);
      }
    }
    return out;
  }



  ROOT::VecOps::RVec<float> XPtoPar_dz(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                                        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
                                        const TLorentzVector& V, // primary vertex
                                        const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e-9; //Reduced speed of light ???

    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        float D0_wrt0 = tracks.at(rp.tracks_begin).D0;
        float Z0_wrt0 = tracks.at(rp.tracks_begin).Z0;
        float phi0_wrt0 = tracks.at(rp.tracks_begin).phi;

        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - V.Vect();

        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);

        double a = - rp.charge * Bz * cSpeed;
        double pt = p.Pt();
        double C = a/(2 * pt);
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
        double D;
        if (pt < 10.0) D = (T - pt) / a;
        else D = (-2 * cross + a * r2) / (T + pt);
        double B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
        if ( TMath::Abs(B) > 1.) B = TMath::Sign(1, B);
        double st = TMath::ASin(B) / C;
        double ct = p(2) / pt;
        double z0;
        double dot = x(0) * p(0) + x(1) * p(1);
        if (dot > 0.0) z0 = x(2) - ct * st;
        else z0 = x(2) + ct * st;

        out.push_back(z0);
      } else {
        out.push_back(-9.);
      }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_phi(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
					const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
					const TLorentzVector& V, // primary vertex
					const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e-9; //Reduced speed of light ???

    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        float D0_wrt0 = tracks.at(rp.tracks_begin).D0;
        float Z0_wrt0 = tracks.at(rp.tracks_begin).Z0;
        float phi0_wrt0 = tracks.at(rp.tracks_begin).phi;

        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - V.Vect();

        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);

        double a = - rp.charge * Bz * cSpeed;
        double pt = p.Pt();
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
        double phi0 = TMath::ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T);

	out.push_back(phi0);

      } else {
        out.push_back(-9.);
      }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_C(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
				       const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
				       const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e3 * 1.0e-15;
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);

        double a = std::copysign(1.0, rp.charge) * Bz * cSpeed;
	double pt = p.Pt();
        double C = a/(2 * pt);

	out.push_back(C);
      } else {
        out.push_back(-9.);
      }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_ct(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
				       const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
				       const float& Bz) {

    const double cSpeed = 2.99792458e8 * 1.0e-9;
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {

      if( rp.tracks_begin < tracks.size()) {

        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
	double pt = p.Pt();

        double ct = p(2) / pt;

	out.push_back(ct);

      } else {
        out.push_back(-9.);
      }
    }
    return out;
  }


ROOT::VecOps::RVec<float>
getRP2TRK_D0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
             ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).D0);
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_D0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					      ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[0]);
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_D0_sig(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					      ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).D0/sqrt(tracks.at(p.tracks_begin).covMatrix[0]));
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_Z0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).Z0);
    else result.push_back(-9.);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_Z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					      ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[9]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_Z0_sig(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					      ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).Z0/sqrt(tracks.at(p.tracks_begin).covMatrix[9]));
    else result.push_back(std::nan(""));
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					   ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).phi);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					       ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[2]);
    else result.push_back(-9);
  }
  return result;
}


ROOT::VecOps::RVec<float>
getRP2TRK_omega(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					     ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).omega);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_omega_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						 ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[5]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_tanLambda(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						 ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).tanLambda);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_tanLambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						     ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[14]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_phi0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						   ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[1]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_omega_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						    ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[3]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						 ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[6]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_d0_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
							ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[10]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi0_omega_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						      ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[4]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi0_z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						   ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[7]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_phi0_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
							  ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[11]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_omega_z0_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
						    ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[8]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_omega_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
							   ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[12]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<float>
getRP2TRK_z0_tanlambda_cov(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
							ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    if (p.tracks_begin<tracks.size())
      result.push_back(tracks.at(p.tracks_begin).covMatrix[13]);
    else result.push_back(-9);
  }
  return result;
}

ROOT::VecOps::RVec<edm4hep::TrackState>
getRP2TRK( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
					ROOT::VecOps::RVec<edm4hep::TrackState> tracks )
{

  ROOT::VecOps::RVec<edm4hep::TrackState> result ;
  result.reserve( in.size() );

  for (auto & p: in) {
    if (p.tracks_begin >= 0 && p.tracks_begin<tracks.size()) {
	result.push_back(tracks.at(p.tracks_begin) ) ;
    }
  }
 return result ;
}

// returns reco indices of tracks
ROOT::VecOps::RVec<int> 
get_recoindTRK( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, 
		ROOT::VecOps::RVec<edm4hep::TrackState> tracks )
{

  ROOT::VecOps::RVec<int> result ;
  
  for (unsigned int ctr=0; ctr<in.size(); ctr++) {
    edm4hep::ReconstructedParticleData p = in.at(ctr);
    if (p.tracks_begin >= 0 && p.tracks_begin<tracks.size()) result.push_back(ctr) ;
  }
 return result ;
}

int getTK_n(ROOT::VecOps::RVec<edm4hep::TrackState> x) {
  int result =  x.size();
  return result;
}

///
ROOT::VecOps::RVec<bool> 
hasTRK( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) {

  ROOT::VecOps::RVec<bool> result ;
  result.reserve( in.size() );
  
  for (auto & p: in) {
    if (p.tracks_begin >= 0 && p.tracks_begin != p.tracks_end) result.push_back(true) ;
    else result.push_back(false);
  }
 return result ;
}

}//end NS ReconstructedParticle2Track

}//end NS FCCAnalyses

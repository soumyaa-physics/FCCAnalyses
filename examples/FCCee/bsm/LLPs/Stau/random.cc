random root 

events->Draw("GenTau_cTau/100:GenTau_theta*180/3.14159", "", "colz");
c1->SaveAs("cTau_vs_theta.png");


/// Kink Finder
ROOT::VecOps::RVec<float>findKink_angle(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced
) {
   // angle between primary and displaced track
    ROOT::VecOps::RVec<float> out;

    // useful for background
    if (displaced.empty()) return out;

    for (size_t i = 0; i < primary.size(); ++i) {
        const auto& tIn = primary[i]; // current primary track
        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) { // get associated reco particle
            if (rp.tracks_begin == i) { rpIn = &rp; break; }
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        // doing the same for displaced tracks
        for (size_t j = 0; j < displaced.size(); ++j) {
            const auto& tOut = displaced[j];
            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                if (rp.tracks_begin == j) { rpOut = &rp; break; }
            }
            if (!rpOut) continue;
            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) > 0.05) continue;
            // define kink angle
            // momentum is defined in vertexing utils - but float
            TVector3 pIn(rpIn->momentum.x, rpIn->momentum.y, rpIn->momentum.z);
            TVector3 pOut(rpOut->momentum.x, rpOut->momentum.y, rpOut->momentum.z);
            float cosTheta = pIn.Dot(pOut) / (pIn.Mag() * pOut.Mag());
            float kinkAngle = std::acos(std::clamp(cosTheta, -1.f, 1.f));
            out.push_back(kinkAngle);
        }
    }
    return out;
  }

  ROOT::VecOps::RVec<int> findKink_size(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced
) {
   // output number of kink candidates
    ROOT::VecOps::RVec<int> out;
    int nKinks = 0;

    // useful for background
    if (displaced.empty()) return out;

    for (size_t i = 0; i < primary.size(); ++i) {
        const auto& tIn = primary[i]; // current primary track
        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) { // get associated reco particle
            if (rp.tracks_begin == i) { rpIn = &rp; break; }
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        
        // doing the same for displaced tracks
        for (size_t j = 0; j < displaced.size(); ++j) {
            const auto& tOut = displaced[j];
            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                if (rp.tracks_begin == j) { rpOut = &rp; break; }
            }
            if (!rpOut) continue;
            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) > 0.05) continue;
            ++nKinks;
          }
          out.push_back(nKinks);
    }
    return out;
  }

    ROOT::VecOps::RVec<float> findKink_x(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced
) {
   // output x coordinate of kink candidates
    ROOT::VecOps::RVec<float> out;

    // useful for background
    if (displaced.empty()) return out;

    for (size_t i = 0; i < primary.size(); ++i) {
        const auto& tIn = primary[i]; // current primary track
        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) { // get associated reco particle
            if (rp.tracks_begin == i) { rpIn = &rp; break; }
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        
        // doing the same for displaced tracks
        for (size_t j = 0; j < displaced.size(); ++j) {
            const auto& tOut = displaced[j];
            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                if (rp.tracks_begin == j) { rpOut = &rp; break; }
            }
            if (!rpOut) continue;
            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) > 0.05) continue;

            float x = tOut.referencePoint.x;
            out.push_back(x);
            }
    }
    return out;
  }

    ROOT::VecOps::RVec<float> findKink_y(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced
) {
   // output y coordinate of kink candidates
    ROOT::VecOps::RVec<float> out;

    // useful for background
    if (displaced.empty()) return out;

    for (size_t i = 0; i < primary.size(); ++i) {
        const auto& tIn = primary[i]; // current primary track
        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) { // get associated reco particle
            if (rp.tracks_begin == i) { rpIn = &rp; break; }
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        
        // doing the same for displaced tracks
        for (size_t j = 0; j < displaced.size(); ++j) {
            const auto& tOut = displaced[j];
            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                if (rp.tracks_begin == j) { rpOut = &rp; break; }
            }
            if (!rpOut) continue;
            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) > 0.05) continue;
            float y = tOut.referencePoint.y;
            out.push_back(y);
        }
    }
    return out;
  }

   ROOT::VecOps::RVec<float> findKink_z(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& primary,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& displaced
) {
   // output z coordinate of kink candidates
    ROOT::VecOps::RVec<float> out;

    // useful for background
    if (displaced.empty()) return out;

    for (size_t i = 0; i < primary.size(); ++i) {
        const auto& tIn = primary[i]; // current primary track
        const edm4hep::ReconstructedParticleData* rpIn = nullptr; 
        for (const auto& rp : reco) { // get associated reco particle
            if (rp.tracks_begin == i) { rpIn = &rp; break; }
        }
        if (!rpIn) continue; // if no reco particle associated to this track, skip
        
        // doing the same for displaced tracks
        for (size_t j = 0; j < displaced.size(); ++j) {
            const auto& tOut = displaced[j];
            const edm4hep::ReconstructedParticleData* rpOut = nullptr;
            for (const auto& rp : reco) {
                if (rp.tracks_begin == j) { rpOut = &rp; break; }
            }
            if (!rpOut) continue;
            // check charge, impact parameters
            if (rpIn->charge != rpOut->charge) continue;
            if (std::abs(tIn.D0 - tOut.D0) > 0.05) continue;
            float z = tOut.referencePoint.z;
            out.push_back(z);
        }
    }
    return out;
  }

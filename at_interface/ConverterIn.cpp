#include "ConverterIn.hpp"

#include "Constants.hpp"
#include "InputContainer.hpp"

#include "AnalysisTree/TaskManager.hpp"
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/Matching.hpp>

using namespace AnalysisTree;

void ConverterIn::FillParticle(const AnalysisTree::BranchChannel& rec_particle) {
  
  std::vector<float> mf(NumberOfFieldPars, 0.f);
  mf.at(0) = rec_particle[kf_tracks_.GetField("cx0")];
  mf.at(1) = rec_particle[kf_tracks_.GetField("cx1")];
  mf.at(2) = rec_particle[kf_tracks_.GetField("cx2")];
  mf.at(3) = rec_particle[kf_tracks_.GetField("cy0")];
  mf.at(4) = rec_particle[kf_tracks_.GetField("cy1")];
  mf.at(5) = rec_particle[kf_tracks_.GetField("cy2")];
  mf.at(6) = rec_particle[kf_tracks_.GetField("cz0")];
  mf.at(7) = rec_particle[kf_tracks_.GetField("cz1")];
  mf.at(8) = rec_particle[kf_tracks_.GetField("cz2")];
  mf.at(9) = rec_particle[kf_tracks_.GetField("z0")];

  auto cov_matrix = is_shine_ ? GetCovMatrixShine(rec_particle) : GetCovMatrixCbm(rec_particle);
  
  std::vector<float> par(kNumberOfTrackPars, 0.f);
  par.at(kX) = rec_particle[kf_tracks_.GetField("x")];
  par.at(kY) = rec_particle[kf_tracks_.GetField("y")];
  par.at(kZ) = rec_particle[kf_tracks_.GetField("z")];
  par.at(kPx) = rec_particle[kf_tracks_.GetField("px")];
  par.at(kPy) = rec_particle[kf_tracks_.GetField("py")];
  par.at(kPz) = rec_particle[kf_tracks_.GetField("pz")];

  const int q = rec_particle[kf_tracks_.GetField("q")];

  int pdg = -999;
  if (pid_mode_ == 0) {
    pdg = rec_particle[kf_tracks_.GetField("q")];
    container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[kf_tracks_.GetField("nhits")]);
  } else if (pid_mode_ == 1) {
    pdg = rec_particle[kf_tracks_.GetField("mc_pdg")];
    container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[kf_tracks_.GetField("nhits")]);
  } else {
    if (rec_particle[kf_tracks_.GetField("prob_p")] == -1)// needs to be done to exclude tracks with no TOF id (tracks with no TOF id have the same pid than negative background)
      return;

    if (pid_mode_ == 2 && pid_purity_.at(0) == 0.5) {
      const int pdg = rec_particle[kf_tracks_.GetField("rec_pdg")] * q;    // Be careful if use electrons and muons, because then the sign of pdg and charge do not match
      container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[kf_tracks_.GetField("nhits")]);

    } else {

      std::vector<float> pdg_prob;
      pdg_prob.push_back(rec_particle[kf_tracks_.GetField("prob_p")]);
      pdg_prob.push_back(rec_particle[kf_tracks_.GetField("prob_pi")]);
      pdg_prob.push_back(rec_particle[kf_tracks_.GetField("prob_K")]);
      pdg_prob.push_back(rec_particle[kf_tracks_.GetField("prob_d")]);
      pdg_prob.push_back(rec_particle[kf_tracks_.GetField("prob_bg")]);      

      if (pid_mode_ == 2) {
        if (*std::max_element(pdg_prob.begin(), pdg_prob.end()) < pid_purity_.at(0))
          return;
        auto it_prob = std::max_element(pdg_prob.begin(), pdg_prob.end());
        int ipid = std::distance(pdg_prob.begin(), it_prob);
        pdg = pid_codes_rec[ipid] * q;
        container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[kf_tracks_.GetField("nhits")]);
      }

      if (pid_mode_ == 3) {
        for (size_t ipid = 0; ipid < pid_codes_rec.size(); ipid++)
          if (pdg_prob[ipid] >= pid_purity_.at(ipid)) {
            pdg = pid_codes_rec[ipid] * q;
            container_.AddTrack(par, cov_matrix, mf, q, pdg, rec_particle.GetId(), rec_particle[kf_tracks_.GetField("nhits")]);
          }
      }
    }
  }
}
void ConverterIn::Init() {
  auto* chain = AnalysisTree::TaskManager::GetInstance()->GetChain();

  if (pid_mode_ > 1) kf_tracks_name_ = "RecTracks";

  rec_event_header_ = chain->GetBranch(rec_event_header_name_);
  sim_event_header_ = chain->GetBranch(sim_event_header_name_);
  kf_tracks_ = chain->GetBranch(kf_tracks_name_);
  sim_tracks_ = chain->GetBranch(sim_tracks_name_);
  kf2sim_tracks_ = chain->GetMatching(kf_tracks_name_, sim_tracks_name_);

// //   const auto& branch_conf_kftr = config_->GetBranchConfig(kf_tracks_name_);
// //   q_field_id_ = branch_conf_kftr.GetFieldId("q");
// //   par_field_id_ = branch_conf_kftr.GetFieldId("x");   // par0
// //   mf_field_id_ = branch_conf_kftr.GetFieldId("cx0");  // magnetic field par0
// //   cov_field_id_ = branch_conf_kftr.GetFieldId("cov1");// cov matrix 0
// //   passcuts_field_id_ = branch_conf_kftr.GetFieldId("pass_cuts");
// //   pdg_field_id_ = branch_conf_kftr.GetFieldId("mc_pdg");
// //   nhits_field_id_ = branch_conf_kftr.GetFieldId("nhits");
// // 
// //   if (pid_mode_ > 1) {
// //     pdg_prob_field_id_ = branch_conf_kftr.GetFieldId("prob_p");
// //     pdg_rec_field_id_ = branch_conf_kftr.GetFieldId("rec_pdg");
// //   }
// // 
// //   const auto& branch_conf_simtr = config_->GetBranchConfig(sim_tracks_name_);
// //   mother_id_field_id_ = branch_conf_simtr.GetFieldId("mother_id");
// //   sim_pdg_field_id_ = branch_conf_simtr.GetFieldId("pdg");

  if (track_cuts_) {
    track_cuts_->Init(*config_);
  }
}

void ConverterIn::Exec() {

  container_ = InputContainer();
  const int n_tracks = kf_tracks_.size();

  container_.SetPV(rec_event_header_[0][rec_event_header_.GetField("vtx_x")],
                   rec_event_header_[0][rec_event_header_.GetField("vtx_y")],
                   rec_event_header_[0][rec_event_header_.GetField("vtx_z")]);

  int n_good_tracks{0};
  container_.Reserve(n_tracks);
  for (int i_track = 0; i_track < n_tracks; ++i_track) {
    const auto& rec_track = kf_tracks_[i_track];
    if (!IsGoodTrack(rec_track)) continue;
    if (!CheckMotherPdgs(rec_track)) continue;
    FillParticle(rec_track);
    n_good_tracks++;
  }
}

std::vector<float> ConverterIn::GetCovMatrixCbm(const AnalysisTree::BranchChannel& particle) const {
  const auto tx = particle[kf_tracks_.GetField("tx")];
  const auto ty = particle[kf_tracks_.GetField("ty")];
  const auto qp = particle[kf_tracks_.GetField("qp")];
  const auto q =  particle[kf_tracks_.GetField("q")];

  //calculate covariance matrix
  const auto t = sqrt(1.f + tx * tx + ty * ty);
  const auto t3 = t * t * t;
  const auto dpxdtx = q / qp * (1.f + ty * ty) / t3;
  const auto dpxdty = -q / qp * tx * ty / t3;
  const auto dpxdqp = -q / (qp * qp) * tx / t;
  const auto dpydtx = -q / qp * tx * ty / t3;
  const auto dpydty = q / qp * (1.f + tx * tx) / t3;
  const auto dpydqp = -q / (qp * qp) * ty / t;
  const auto dpzdtx = -q / qp * tx / t3;
  const auto dpzdty = -q / qp * ty / t3;
  const auto dpzdqp = -q / (qp * qp * t3);

  const float F[kNumberOfTrackPars][5] = {{1.f, 0.f, 0.f, 0.f, 0.f},
                                          {0.f, 1.f, 0.f, 0.f, 0.f},
                                          {0.f, 0.f, 0.f, 0.f, 0.f},
                                          {0.f, 0.f, dpxdtx, dpxdty, dpxdqp},
                                          {0.f, 0.f, dpydtx, dpydty, dpydqp},
                                          {0.f, 0.f, dpzdtx, dpzdty, dpzdqp}};

  float VFT[5][kNumberOfTrackPars];
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < kNumberOfTrackPars; j++) {
      VFT[i][j] = 0;
      for (int k = 0; k < 5; k++) {
        VFT[i][j] += particle[kf_tracks_.GetField(("cov" + std::to_string(1 + std::min(i, k) + std::max(i, k) * (std::max(i, k) + 1) / 2)).c_str())] * F[j][k];
      }
    }
  }

  std::vector<float> cov(21, 0);
  for (int i = 0, l = 0; i < kNumberOfTrackPars; i++) {
    for (int j = 0; j <= i; j++, l++) {
      cov.at(l) = 0;
      for (int k = 0; k < 5; k++) {
        cov.at(l) += F[i][k] * VFT[k][j];
      }
    }
  }
  return cov;
}

std::vector<float> ConverterIn::GetCovMatrixShine(const AnalysisTree::BranchChannel& particle) const {
  std::vector<float> cov(21, 0.);
  for (int iCov = 0; iCov < 21; ++iCov) {
    cov[iCov] = particle[kf_tracks_.GetField(("cov" + std::to_string(iCov+1)).c_str())];
  }
  return cov;
}

bool ConverterIn::IsGoodTrack(const AnalysisTree::BranchChannel& rec_track) const {
  return track_cuts_ ? track_cuts_->Apply(rec_track) : true;
}

bool ConverterIn::CheckMotherPdgs(const AnalysisTree::BranchChannel& rec_track) const {
  if (mother_pdgs_to_be_considered_.size() == 0)
    return true;

//   if (!sim_tracks_ || !kf2sim_tracks_) {
//     std::cout << "No MC info available!\n";
//     assert(false);
//   }

  const int sim_id = kf2sim_tracks_->GetMatch(rec_track.GetId());
  if (sim_id < 0)
    return false;

  const AnalysisTree::BranchChannel& sim_track = sim_tracks_[sim_id];
  const int mother_id = sim_track[sim_tracks_.GetField("mother_id")];
  if (mother_id < 0)
    return false;

  const int mother_pdg = sim_tracks_[mother_id][sim_tracks_.GetField("pid")];

  bool ok = false;

  for (auto& good_mother_pdgs : mother_pdgs_to_be_considered_)
    if (mother_pdg == good_mother_pdgs) {
      ok = true;
      break;
    }

  return ok;
}

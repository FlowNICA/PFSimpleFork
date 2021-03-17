#include <KFSimple/Constants.h>
#include <Interface/InputContainer.h>
#include "ConverterIn.h"

void ConverterIn::FillParticle(const AnalysisTree::Track& rec_particle, InputContainer& input_info) const
{
  std::vector <float> mf(kNumberOfFieldPars, 0.f);
  for(int iF=0; iF<kNumberOfFieldPars; iF++){
    mf[iF] = rec_particle.GetField<float>(mf_field_id_+iF);
  }
  auto cov_matrix = is_shine_ ? GetCovMatrixShine(rec_particle) : GetCovMatrixCbm(rec_particle);
  std::vector <float> par(kNumberOfTrackPars,0.f);

  par[kX] = rec_particle.GetField<float>(par_field_id_);
  par[kY] = rec_particle.GetField<float>(par_field_id_ + 1);
  par[kZ] = rec_particle.GetField<float>(par_field_id_ + 2);
  par[kPx] = rec_particle.GetPx();
  par[kPy] = rec_particle.GetPy();
  par[kPz] = rec_particle.GetPz();

  const int pdg = rec_particle.GetField<int>(pdg_field_id_); //TODO
//  const int pdg = rec_particle.GetPid();

  input_info.AddTrack(par, cov_matrix, mf, rec_particle.GetField<int>(q_field_id_), pdg, rec_particle.GetField<float>(chi2_field_id_), rec_particle.GetField<int>(ndf_field_id_), rec_particle.GetId(), rec_particle.GetField<int>(nhits_field_id_));
}


void ConverterIn::Init(std::map<std::string, void*>& branches) {
  rec_event_header_ = (AnalysisTree::EventHeader*) branches.find(in_branches_[kRecEventHeader])->second;
  sim_event_header_ = (AnalysisTree::EventHeader*) branches.find(in_branches_[kSimEventHeader])->second;
  kf_tracks_ = (AnalysisTree::TrackDetector*) branches.find(in_branches_[kKfpfTracks])->second;
  sim_tracks_ = (AnalysisTree::Particles*) branches.find(in_branches_[kSimTracks])->second;
  kf2sim_tracks_ = (AnalysisTree::Matching*) branches.find(
    config_->GetMatchName(in_branches_[kKfpfTracks], in_branches_[kSimTracks]))->second;

  auto branch_conf_kftr = config_->GetBranchConfig(in_branches_[kKfpfTracks]);

  q_field_id_ = branch_conf_kftr.GetFieldId("q");
  par_field_id_ = branch_conf_kftr.GetFieldId("x");   // par0
  mf_field_id_ = branch_conf_kftr.GetFieldId("cx0");  // magnetic field par0
  cov_field_id_ = branch_conf_kftr.GetFieldId("cov1"); // cov matrix 0
  passcuts_field_id_ = branch_conf_kftr.GetFieldId("pass_cuts");
  pdg_field_id_ = branch_conf_kftr.GetFieldId("mc_pdg");
  nhits_field_id_ = branch_conf_kftr.GetFieldId("nhits");
  chi2_field_id_ = branch_conf_kftr.GetFieldId("chi2");
  ndf_field_id_ = branch_conf_kftr.GetFieldId("ndf");

  auto branch_conf_simtr = config_->GetBranchConfig(in_branches_[kSimTracks]);
  mother_id_field_id_ = branch_conf_simtr.GetFieldId("mother_id");
  sim_pdg_field_id_ = branch_conf_simtr.GetFieldId("pdg");

  if (track_cuts_){
    track_cuts_->Init(*config_);
  }
}

InputContainer ConverterIn::CreateInputContainer() const {

  InputContainer input_container;
  const int n_tracks = kf_tracks_->GetNumberOfChannels();
//  std::cout << " Ntracks = " << n_tracks << std::endl;
  input_container.SetCuts(cuts_);
  input_container.SetPV(rec_event_header_->GetVertexX(), rec_event_header_->GetVertexY(), rec_event_header_->GetVertexZ());

  int n_good_tracks{0};
  input_container.Reserve(n_tracks);
  for(int i_track=0; i_track<n_tracks; ++i_track) {
    const auto& rec_track = kf_tracks_->GetChannel(i_track);
    if(!IsGoodTrack(rec_track)) continue;
    FillParticle(rec_track, input_container);
    n_good_tracks++;
  }
//  std::cout << "Good tracks = " << n_good_tracks << "\n" << std::endl;
  return input_container;
}

std::vector<float> ConverterIn::GetCovMatrixCbm(const AnalysisTree::Track& particle) const {
  const auto tx = particle.GetField<float>(par_field_id_ + 3);
  const auto ty = particle.GetField<float>(par_field_id_ + 4);
  const auto qp = particle.GetField<float>(par_field_id_ + 5);
  const auto q = particle.GetField<int>(q_field_id_);

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
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < kNumberOfTrackPars; j++) {
      VFT[i][j] = 0;
      for(int k = 0; k < 5; k++) {
        VFT[i][j] += particle.GetField<float>(cov_field_id_ + std::min(i, k) + std::max(i, k) * (std::max(i, k) + 1) / 2) * F[j][k];   //parameters->GetCovariance(i,k) * F[j][k];
//        if(k <= i)
//          VFT[i][j] += particle.GetField<float>(cov_field_id_ + k + i * (i + 1) / 2) * F[j][k];   //parameters->GetCovariance(i,k) * F[j][k];
//        else
//          VFT[i][j] += particle.GetField<float>(cov_field_id_ + i + k * (k + 1) / 2) * F[j][k];   //parameters->GetCovariance(i,k) * F[j][k];
      }
    }
  }

  std::vector<float> cov(21, 0);
  for(int i=0, l=0; i<kNumberOfTrackPars; i++){
    for(int j=0; j<=i; j++, l++) {
      cov[l] = 0;
      for(int k=0; k<5; k++) {
        cov[l] += F[i][k] * VFT[k][j];
      }
    }
  }
  return cov;
}

std::vector<float> ConverterIn::GetCovMatrixShine(const AnalysisTree::Track& particle) const
{
  std::vector<float> cov(21, 0.);

  for (int iCov=0; iCov<21; ++iCov)
    cov[iCov] = particle.GetField<float>(cov_field_id_+iCov);

  return cov;
}

bool ConverterIn::IsGoodTrack(const AnalysisTree::Track& rec_track) const
{
  if(!track_cuts_){
    return true;
  }
  return track_cuts_->Apply(rec_track);

//  return true;
  bool is_good{false};

  const int sim_id = kf2sim_tracks_->GetMatch(rec_track.GetId());
//    std::cout<< "sim  "  << sim_id << std::endl;
  if (sim_id >= 0 && sim_id < sim_tracks_->GetNumberOfChannels()) {
    const AnalysisTree::Track& sim_track = sim_tracks_->GetChannel(sim_id);
    const int mother_id = sim_track.GetField<int>(mother_id_field_id_);
//      std::cout << "mother " << mother_id << std::endl;

    if (mother_id >= 0 && mother_id < sim_tracks_->GetNumberOfChannels() )
    {
      const AnalysisTree::Track& mother_track = sim_tracks_->GetChannel(mother_id);
      const int mother_pdg = mother_track.GetField<int>(sim_pdg_field_id_);
      std::cout << "mother pdg " << mother_pdg << std::endl;

      if(mother_pdg == 3122)
        is_good = true;
    }
  }
  return is_good;
}

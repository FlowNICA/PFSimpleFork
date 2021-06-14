#include "ConverterOut.hpp"

#include "PFSimpleTask.hpp"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <algorithm>

void ConverterOut::CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) const {

  particle.SetMomentum(kf_particle.GetPx(), kf_particle.GetPy(), kf_particle.GetPz());
  particle.SetMass(kf_particle.GetMass());
  particle.SetPid(kf_particle.GetPdg());

  particle.SetField(kf_particle.GetX(), x_field_id_);
  particle.SetField(kf_particle.GetY(), x_field_id_ + 1);
  particle.SetField(kf_particle.GetZ(), x_field_id_ + 2);
  particle.SetField(kf_particle.GetXError(), x_field_id_ + 3);
  particle.SetField(kf_particle.GetYError(), x_field_id_ + 4);
  particle.SetField(kf_particle.GetZError(), x_field_id_ + 5);

  particle.SetField(kf_particle.GetPtError(), pt_err_field_id_);
  particle.SetField(kf_particle.GetPhiError(), pt_err_field_id_ + 1);
  particle.SetField(kf_particle.GetEtaError(), pt_err_field_id_ + 2);
  particle.SetField(kf_particle.GetMassError(), pt_err_field_id_ + 3);
  
  const int Ndaughters = kf_particle.GetDaughterIds().size();

  for (int i = 0; i < Ndaughters; ++i) {
    particle.SetField(kf_particle.GetChi2Prim(i), chi2prim_field_id_ + i);
    particle.SetField(kf_particle.GetCos(i), cosine_field_id_ + i);
    particle.SetField(kf_particle.GetDaughterIds().at(i), daughter_id_field_id_+i);
    particle.SetField(kf_particle.GetInvMassDiscr(i), invmass_discr_field_id_+i);
  }
  
  if(Ndaughters < 3) {
    particle.SetField(-999.f, chi2prim_field_id_ + 2);
    particle.SetField(-9.f, cosine_field_id_ + 2);
    particle.SetField(-999, daughter_id_field_id_ + 2);
    particle.SetField(-999.f, invmass_discr_field_id_ + 2);
    particle.SetField(-999.f, distance_field_id_ + 1);
  }

  particle.SetField(kf_particle.GetDistance(), distance_field_id_);

  particle.SetField(kf_particle.GetChi2Geo(), chi2geo_field_id_);
  particle.SetField(kf_particle.GetL(), chi2geo_field_id_ + 1);
  particle.SetField(kf_particle.GetLdL(), chi2geo_field_id_ + 2);
  particle.SetField(kf_particle.GetChi2Topo(), chi2geo_field_id_ + 3);
  particle.SetField(kf_particle.GetCosineTopo(), chi2geo_field_id_ + 4);
}

void ConverterOut::Exec() {
  
  candidates_ = pfsimple_task_->GetSimpleFinder()->GetCandidates();
  
  lambda_reco_->ClearChannels();
  lambda_sim_->ClearChannels();
  lambda_reco2sim_->Clear();

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  events_->SetField(float(sim_events_->GetField<float>(b_field_id_)), 0);//TODO

  const auto& br_conf = out_config->GetBranchConfig(lambda_reco_->GetId());

  for (const auto& candidate : candidates_) {
    auto& lambdarec = lambda_reco_->AddChannel(br_conf);
    CopyParticle(candidate, lambdarec);
  }
  
  if(mc_particles_ != nullptr)
    MatchWithMc();
}

void ConverterOut::Init() {
  
  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();

  sim_events_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(sim_events_name_));
  mc_particles_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch(mc_particles_name_));
  rec_tracks_ = ANALYSISTREE_UTILS_GET<AnalysisTree::TrackDetector*>(chain->GetPointerToBranch(rec_tracks_name_));
  rec_to_mc_ = chain->GetMatchPointers().find(config_->GetMatchName(rec_tracks_name_, mc_particles_name_))->second;

  std::string out_branch_event = "Events";
  std::string out_branch = std::string("Candidates");
  std::string out_branch_sim = std::string("Simulated");
  std::string out_branch_reco2sim = out_branch + "2" + out_branch_sim;

  AnalysisTree::BranchConfig EventBranch(out_branch_event, AnalysisTree::DetType::kEventHeader);
  EventBranch.AddField<float>("b");

  AnalysisTree::BranchConfig out_particles(out_branch, AnalysisTree::DetType::kParticle);

  out_particles.AddFields<float>({"x", "y", "z", "x_error", "y_error", "z_error"});
  out_particles.AddFields<float>({"pT_err", "phi_err", "eta_err", "mass_err"});

  out_particles.AddFields<int>({"daughter1_id", "daughter2_id", "daughter3_id"});
  out_particles.AddFields<float>({"chi2_prim_first", "chi2_prim_second", "chi2_prim_third"});
  out_particles.AddFields<float>({"distance", "distance_third"});
  out_particles.AddFields<float>({"cosine_first", "cosine_second", "cosine_third"});

  out_particles.AddFields<float>({"chi2_geo", "l", "l_over_dl", "chi2_topo", "cosine_topo"});
  out_particles.AddFields<float>({"invmass_discr_first", "invmass_discr_second", "invmass_discr_third"});

  if (mc_particles_) {
    out_particles.AddField<int>("generation");
  }
  AnalysisTree::BranchConfig LambdaSimBranch(out_branch_sim, AnalysisTree::DetType::kParticle);
  LambdaSimBranch.AddField<int>("mother_id", "particle mother's id in SimParticles branch");   // particle mother's id in SimParticles branch

  man->AddBranch(out_branch_event, events_, EventBranch);
  man->AddBranch(out_branch, lambda_reco_, out_particles);
  man->AddBranch(out_branch_sim, lambda_sim_, LambdaSimBranch);
  man->AddMatching(out_branch, out_branch_sim, lambda_reco2sim_);

  events_->Init(EventBranch);
  InitIndexes();
}

int ConverterOut::GetMothersSimId(AnalysisTree::Particle& lambdarec)
{
  /** Returnes mother's id in SimParticles branch if:
   ** 1) all reconstructed daughters have matching with simulated particles
   ** 2) all of them belong to the same mother
   ** 3) mother has PDG code as expected according to decay reconstruction hypothesis
   **/ 
  int Ndaughters;
  if(lambdarec.GetField<int>(daughter_id_field_id_+2) == -999)
    Ndaughters = 2;
  else
    Ndaughters = 3;
  
  std::vector<int> daughter_sim_id;
  
  if(lambdarec.GetPid() == 3122) {            // TODO make flag (or smth else or universalize) whether decay is cascade or not
    for(int i=0; i<Ndaughters; i++)
      daughter_sim_id.push_back(rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_+i)));
    
    if(*std::min_element(daughter_sim_id.begin(), daughter_sim_id.end())<0) // at least one daughter has no matching with mc
      return -1; 
    
    std::vector<int> mother_sim_id;
    for(int i=0; i<Ndaughters; i++)
      mother_sim_id.push_back(mc_particles_->GetChannel(daughter_sim_id.at(i)).GetField<int>(mother_id_field_id_));
    
    if(*std::min_element(mother_sim_id.begin(), mother_sim_id.end()) != *std::max_element(mother_sim_id.begin(), mother_sim_id.end())) // daughters belong to not the same mother
      return -1;
    
    if(mother_sim_id.at(0)<0) // mother has negative id
      return -1;
    
    if(mc_particles_->GetChannel(mother_sim_id.at(0)).GetPid() != lambdarec.GetPid()) // mother has not PDG which was supposed
      return -1;
      
    return mother_sim_id.at(0);
  }
  else if(lambdarec.GetPid() == 3312) {      // For Xi- the first daughter is pion, the second - Lambda
    daughter_sim_id.push_back(rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_+0)));
    daughter_sim_id.push_back(lambda_reco2sim_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_+1)));
    
//     std::cout << "lambdarec.GetField<int>(daughter_id_field_id_+1) = " << lambdarec.GetField<int>(daughter_id_field_id_+1) << "\tlambda_reco2sim_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_+1)) = " << lambda_reco2sim_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_+1)) << "\n";
    
    if(*std::min_element(daughter_sim_id.begin(), daughter_sim_id.end())<0) // at least one daughter has no matching with mc
      return -1;
    
    std::vector<int> mother_sim_id;
    mother_sim_id.push_back(mc_particles_->GetChannel(daughter_sim_id.at(0)).GetField<int>(mother_id_field_id_));
    mother_sim_id.push_back(lambda_sim_->GetChannel(daughter_sim_id.at(1)).GetField<int>(mother_id_field_id_w_));
    
    if(*std::min_element(mother_sim_id.begin(), mother_sim_id.end()) != *std::max_element(mother_sim_id.begin(), mother_sim_id.end())) // daughters belong to not the same mother
      return -1;
    
    if(mother_sim_id.at(0)<0) // mother has negative id
      return -1;
    
    if(mc_particles_->GetChannel(mother_sim_id.at(0)).GetPid() != lambdarec.GetPid()) // mother has not PDG which was supposed
      return -1;
        
    return mother_sim_id.at(0);
  }
}

int ConverterOut::DetermineGeneration(int mother_sim_id)
{
  int generation = 0;
  int older_id = mother_sim_id;
  while(older_id>=0)
  {
    const auto& simtrackolder = mc_particles_->GetChannel(older_id);
    older_id = simtrackolder.GetField<int>(mother_id_field_id_);
    generation++;
  }
  
  return generation;
}

void ConverterOut::MatchWithMc() {
  
  for (auto& lambdarec : *lambda_reco_) {
    auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();
    
    int mother_id = GetMothersSimId(lambdarec);
    int generation = DetermineGeneration(mother_id);
    lambdarec.SetField(generation, generation_field_id_);
    
//     if(generation > 1 && lambdarec.GetPid() == 3122)
//     {
//       const int mpdg = mc_particles_->GetChannel(mc_particles_->GetChannel(mother_id).GetField<int>(mother_id_field_id_)).GetPid();
//       if (mpdg == 3312)
//         std::cout << mpdg << "\n";
//     }
    
    if (generation<1) continue;
        
    const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

    auto& lambdasim = lambda_sim_->AddChannel(out_config->GetBranchConfig(lambda_sim_->GetId()));

    lambdasim.SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
    lambdasim.SetMass(simtrackmother.GetMass());
    lambdasim.SetPid(simtrackmother.GetPid());
    lambdasim.SetField(simtrackmother.GetField<int>(mother_id_field_id_), mother_id_field_id_w_);
//     std::cout << "lambdarec.GetId() = " << lambdarec.GetId() << "\tlambdasim.GetId() = " << lambdasim.GetId() << "\n";
    lambda_reco2sim_->AddMatch(lambdarec.GetId(), lambdasim.GetId());    
  }
}

void ConverterOut::InitIndexes() {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  const auto& out_branch = out_config->GetBranchConfig(lambda_reco_->GetId());

  auto branch_conf_sim_event = config_->GetBranchConfig(sim_events_name_);
  b_field_id_ = branch_conf_sim_event.GetFieldId("b");

  x_field_id_ = out_branch.GetFieldId("x");
  daughter_id_field_id_ = out_branch.GetFieldId("daughter1_id");
  pt_err_field_id_ = out_branch.GetFieldId("pT_err");

  if (mc_particles_) {
    auto branch_conf_sim = config_->GetBranchConfig(mc_particles_name_);
    auto out_branch_sim = out_config->GetBranchConfig(lambda_sim_->GetId());
    mother_id_field_id_ = branch_conf_sim.GetFieldId("mother_id");
    generation_field_id_ = out_branch.GetFieldId("generation");
    mother_id_field_id_w_ = out_branch_sim.GetFieldId("mother_id");
  }

  chi2prim_field_id_ = out_branch.GetFieldId("chi2_prim_first");
  distance_field_id_ = out_branch.GetFieldId("distance");
  cosine_field_id_ = out_branch.GetFieldId("cosine_first");

  chi2geo_field_id_ = out_branch.GetFieldId("chi2_geo");
  invmass_discr_field_id_ = out_branch.GetFieldId("invmass_discr_first");
}

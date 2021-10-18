#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

void PFSimpleTask::Exec() {
  pf_simple_ = new SimpleFinderNew();

  pf_simple_->SetDecays(decays_);
  pf_simple_->SetPidMode(pid_mode_);
  if (pid_mode_ == 1) pf_simple_->SetPidPurity(pid_purity_);
  pf_simple_->Init(in_task_->GetInputContainer());
  pf_simple_->FindParticles();
}

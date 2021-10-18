#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_

#include "Decay.hpp"
#include "SimpleFinderNew.hpp"

#include "AnalysisTree/Task.hpp"

class ConverterIn;

class PFSimpleTask : public AnalysisTree::Task {

 public:
  PFSimpleTask() = default;
  ~PFSimpleTask() override = default;

  void Init() override{};
  void Exec() override;
  void Finish() override{};
  void SetInTask(ConverterIn* in_task) { in_task_ = in_task; }
  SimpleFinderNew* GetSimpleFinder() { return pf_simple_; }
  void SetDecays(const std::vector<Decay>& decays) { decays_ = decays; }
  void AddDecay(const Decay& decay) { decays_.emplace_back(decay); }
  void SetPidMode(const int pid_mode) { pid_mode_ = pid_mode; }
  void SetPidPurity(const float pid_purity) { pid_purity_ = pid_purity; }
  
 protected:
  ConverterIn* in_task_{nullptr};
  SimpleFinderNew* pf_simple_{nullptr};
  std::vector<Decay> decays_{};
  int pid_mode_{0};
  float pid_purity_{0.5};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_

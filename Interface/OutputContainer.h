/**
 ** @class OutputContainer
 ** @brief Container with output information about reconstructed particles and geometrical decay parameters (quantities to be cut in order to select particles)
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov
 **
 ** Each particle candidate is characterized with set of geometrical decay parameters. Depending on the
 ** value of each parameter the candidate is saved or rejected.
 ** In order to save the reconstructed particle, the KFParticle object is used. It contains all
 ** information about the particle (mass, momentum etc), and access to this information is
 ** possible via KFParticle methods.
 **/ 

#ifndef OutputContainer_H
#define OutputContainer_H

#include "KFParticle.h"

class OutputContainer
{
  
 public:
  
  OutputContainer() = default;
  virtual ~OutputContainer() = default;  
  
  void SetChi2Prim(float chi2prim){ chi2_prim_ = chi2prim; };
  void SetId(int id){ id_ = id; };
  void SetCharge(int charge){ charge_ = charge; };
  
  float GetChi2Prim() const { return chi2_prim_; };
  int   GetId() const { return id_; };
  int   GetCharge() const { return charge_; };

 protected:
   
  float chi2_prim_{-999.};
  int id_{-999};
  int charge_{-999};
  
};

#endif // OutputContainer_H
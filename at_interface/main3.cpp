#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"
#include "ConverterOutTree.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

int main(int argc, char** argv) {

  if (argc < 2) {
    std::cout << "Wrong number of arguments! Please use:\n  ./main filelist.txt\n";
    return EXIT_FAILURE;
  }

  int pid_mode = 1;
  float pid_purity = 0.5;
  
  const std::string& filename = argv[1];
  
  std::vector<Daughter> daughters = {2212, -211, 1000010020};
  for (int idaughter = 0; idaughter < daughters.size(); idaughter++) {
    daughters.at(idaughter).CancelCuts();
    daughters.at(idaughter).SetCutChi2Prim(18.42);
    //daughters.at(idaughter).SetCutCos(0.99825);
  }
  if(pid_mode == 2) {
    daughters.at(0).SetPidPurity(0.5);
    daughters.at(1).SetPidPurity(0.5);
    daughters.at(2).SetPidPurity(0.3);
  }
  Mother mother(3004);
  mother.CancelCuts();
  mother.SetCutDistance(1.0);
  mother.SetCutChi2GeoSM({3.0});
  //mother.SetCutChi2TopoSM({29});
  mother.SetCutDistanceToSV(1.0);
  //mother.SetCutChi2Geo(3.0);
  //mother.SetCutChi2Topo(29);
  //mother.SetCutLdL(3);
  Decay decay("H3L", mother, {daughters});

  std::string outname;                                                                                   
  if (argc > 2) outname = argv[2];                                                                             
  else outname = "PFSimpleOutput";                                                                         
  std::string outfilename;
  if (pid_mode == 0) outfilename = outname + "_pid_mc.root";               
  else outfilename = outname + "_pid_rec.root";
  
  auto* man = AnalysisTree::TaskManager::GetInstance();
  //man->SetOutputName(outfilename, "pTree");

  std::string tracks_name;
  if (pid_mode == 0) tracks_name = "VtxTracks";
  else tracks_name = "RecTracks";
  
  auto* in_converter = new ConverterIn();
  //in_converter->SetTrackCuts(new AnalysisTree::Cuts("Cut to reproduce KFPF", {AnalysisTree::EqualsCut("VtxTracks.pass_cuts", 1)}));
  //in_converter->SetTrackCuts(new Cuts("Cut to reproduce KFPF", {EqualsCut("VtxTracks.pass_cuts", 1)}));
  in_converter->SetIsShine(false);//TODO maybe change name
  in_converter->SetPidMode(pid_mode);
  in_converter->SetInputBranchNames({"SimParticles", tracks_name, "SimEventHeader", "RecEventHeader"});
  
  auto* pf_task = new PFSimpleTask();
  pf_task->SetInTask(in_converter);
  pf_task->SetDecays({decay});
  pf_task->SetPidMode(pid_mode);
  if (pid_mode == 1) pf_task->SetPidPurity(pid_purity);
  
  auto* out_converter = new ConverterOutTree();
  out_converter->SetPFSimpleTask(pf_task);
  out_converter->SetPidMode(pid_mode);
  out_converter->SetInputBranchNames({"SimParticles", tracks_name, "SimEventHeader", "RecEventHeader"});
  out_converter->SetDecay(decay);
  out_converter->SetOutFilename(outfilename);
  
  man->AddTask(in_converter);
  man->AddTask(pf_task);
  man->AddTask(out_converter);

  if (pid_mode == 0) 
    man->Init({filename}, {"rTree"});
  else
    man->Init({filename}, {"pTree"});
  
  man->Run(-1);// -1 = all events
  man->Finish();
  man->ClearTasks();

  return EXIT_SUCCESS;
}

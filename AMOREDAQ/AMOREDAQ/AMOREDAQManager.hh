#pragma once

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "AMOREAlgs/AMOREChunkFIFO.hh"
#include "AMORESystem/AMORETCB.hh"
#include "DAQ/CupDAQManager.hh"
#include "DAQConfig/AbsConf.hh"
#include "DAQConfig/AbsConfList.hh"
#include "AMOREHDF5/AMOREEDM.hh"
#include "OnlConsts/onlconsts.hh"
#include "DAQUtils/ConcurrentDeque.hh"

class AMOREDAQManager : public CupDAQManager {
public:
  AMOREDAQManager();
  ~AMOREDAQManager() override;

  bool AddADC(AbsConfList * conflist) override;
  bool PrepareDAQ() override;

  virtual void Run();

private:
  void RC_AMORETCB();
  void RC_AMOREDAQ();

  bool ReadConfig();

  template <typename T>
  void FillConfigArray(YAML::Node node, int nch, std::function<void(int, T)> setter,
                       bool inc = false);

  void ReadConfigTCB(YAML::Node ymlnode);
  void ReadConfigADC(YAML::Node ymlnode);

  void TF_ReadData_AMORE();
  void TF_StreamData();
  void TF_SWTrigger(int n);
  void TF_WriteEvent_AMORE();

private:
  AMORETCB & fTCB;

  std::vector<std::unique_ptr<AMOREChunkFIFO>> fFIFOs;
  ConcurrentDeque<Crystal_t> fTriggeredCrystals;

  PROCSTATE fStreamStatus;

  unsigned long fTimeDelta;

  ClassDef(AMOREDAQManager, 0)
};

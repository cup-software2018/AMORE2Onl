#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "AMOREHDF5/AMOREEDM.hh"
#include "HDF5Utils/AbsH5Event.hh"

#include "hdf5.h"

class H5AMOREEvent : public AbsH5Event {
public:
  H5AMOREEvent();
  ~H5AMOREEvent() override;

  void Open() override;
  void Close() override;

  // Append a global event containing multiple Crystal hits
  herr_t AppendEvent(const EventInfo_t & info, const std::vector<Crystal_t> & data);
  herr_t ReadEvent(int n) override;

  void SetNDP(int ndp);
  int GetNDP();

  // Accessor for the currently read event data
  Crystal_t * GetData();

protected:
  herr_t FlushBuffer() override;

private:
  // HDF5 Types and Datasets
  hid_t fChType{H5I_INVALID_HID};
  hid_t fDsetIndex{H5I_INVALID_HID};
  hid_t fDsetChs{H5I_INVALID_HID};
  hid_t fDsetPhonon{H5I_INVALID_HID}; // Separated dataset for phonon
  hid_t fDsetPhoton{H5I_INVALID_HID}; // Separated dataset for photon

  std::uint64_t fTotalCrystals{0};
  int fNDP{0};

  // Write & Read internal buffers for Crystals
  std::vector<CrystalHeader_t> fChBuf;
  std::vector<std::uint16_t> fPhononBuf;
  std::vector<std::uint16_t> fPhotonBuf;

  // Pre-allocated vector buffer for zero-overhead user access
  std::vector<Crystal_t> fDataBuf;

  // Trackers and Cached Dataspaces
  hid_t fCurrentReadFid{H5I_INVALID_HID};
  hid_t fFileSpaceInfo{H5I_INVALID_HID};
  hid_t fFileSpaceIndex{H5I_INVALID_HID};
  hid_t fFileSpaceChs{H5I_INVALID_HID};
  hid_t fFileSpacePhonon{H5I_INVALID_HID};
  hid_t fFileSpacePhoton{H5I_INVALID_HID};

  // ==========================================
  // Full Prefetching Buffers (Memory Safe Window)
  // ==========================================
  std::vector<EventInfo_t> fReadBufInfo;
  std::vector<std::uint64_t> fReadBufIndex;
  std::vector<CrystalHeader_t> fPrefetchChs;
  std::vector<std::uint16_t> fPrefetchPhonon;
  std::vector<std::uint16_t> fPrefetchPhoton;

  int fReadBufStart{-1};
  int fReadBufSize{0};
  std::uint64_t fPrefetchChStart{0};

  ClassDef(H5AMOREEvent, 0)
};

// === inline definitions ===

inline void H5AMOREEvent::SetNDP(int ndp) { fNDP = ndp; }

inline Crystal_t * H5AMOREEvent::GetData() { return fDataBuf.data(); }
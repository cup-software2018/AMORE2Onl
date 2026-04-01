#pragma once

#include <cstdint>
#include <vector>

#include "AMOREHDF5/AMOREEDM.hh"
#include "HDF5Utils/AbsH5Hit.hh"

#include "hdf5.h"

class H5AMOREHit : public AbsH5Hit {
public:
  H5AMOREHit();
  ~H5AMOREHit() override;

  void Open() override;
  void Close() override;

  // Append a single self-triggered Crystal hit
  herr_t AppendHit(const Crystal_t & hit);

  // Read the n-th hit from the stream
  herr_t ReadHit(std::uint64_t n) override;

  void SetNDP(int ndp);
  int GetNDP();

  // Accessor for currently read hit data
  const Crystal_t & GetHit() const;

protected:
  herr_t FlushBuffer() override;

private:
  int fNDP{0};
  hid_t fDsetPhonon{H5I_INVALID_HID};
  hid_t fDsetPhoton{H5I_INVALID_HID};

  // Write buffers for Crystal hits
  std::vector<CrystalHeader_t> fChBuf;
  std::vector<std::uint16_t> fPhononBuf;
  std::vector<std::uint16_t> fPhotonBuf;

  // Data holder for currently read single hit
  Crystal_t fCurrentHit{};

  // Trackers and Cached Dataspaces for reading
  hid_t fCurrentReadFid{H5I_INVALID_HID};
  hid_t fFileSpaceChs{H5I_INVALID_HID};
  hid_t fFileSpacePhonon{H5I_INVALID_HID};
  hid_t fFileSpacePhoton{H5I_INVALID_HID};

  // ==========================================
  // Prefetching Buffers (Hit-based sliding window)
  // ==========================================
  std::vector<CrystalHeader_t> fPrefetchChs;
  std::vector<std::uint16_t> fPrefetchPhonon;
  std::vector<std::uint16_t> fPrefetchPhoton;

  std::uint64_t fReadBufStart{0};
  std::uint64_t fReadBufSize{0};

  ClassDef(H5AMOREHit, 0)
};

// === inline definitions ===

inline void H5AMOREHit::SetNDP(int ndp) { fNDP = ndp; }

inline const Crystal_t & H5AMOREHit::GetHit() const { return fCurrentHit; }
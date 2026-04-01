#include <algorithm>
#include <cstring>

#include "AMOREHDF5/H5AMOREEvent.hh"

ClassImp(H5AMOREEvent)

H5AMOREEvent::H5AMOREEvent()
  : AbsH5Event()
{
}

H5AMOREEvent::~H5AMOREEvent() {}

void H5AMOREEvent::Open()
{
  fEvtType = EventInfo_t::BuildType();
  fChType = CrystalHeader_t::BuildType();

  InitSubRun();

  if (!fWriteTag) { return; }

  if (fFile < 0) {
    Error("Open", "invalid file id (fFile = %d). SetFileId must be called before Open().",
          static_cast<int>(fFile));
    return;
  }

  if (fNDP <= 0 || fNDP > kH5AMORENDPMAX) {
    Error("Open", "Invalid NDP: %d (max %d)", fNDP, kH5AMORENDPMAX);
    return;
  }

  // Create /events group
  {
    htri_t gexists = H5Lexists(fFile, "/events", H5P_DEFAULT);
    if (gexists < 0) {
      Error("Open", "H5Lexists(/events) failed");
      return;
    }
    if (gexists == 0) {
      hid_t grp_events = H5Gcreate2(fFile, "/events", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (grp_events < 0) {
        Error("Open", "Failed to create group /events");
        return;
      }
      H5Gclose(grp_events);
    }
  }

  // Create extendable 1D datasets for metadata
  hsize_t dims[1] = {0};
  hsize_t maxdims[1] = {H5S_UNLIMITED};
  hid_t space1d = H5Screate_simple(1, dims, maxdims);

  hid_t dcpl1d = H5Pcreate(H5P_DATASET_CREATE);
  hsize_t chunk1d = (fBufCap > 0) ? static_cast<hsize_t>(fBufCap) : 1024;
  H5Pset_chunk(dcpl1d, 1, &chunk1d);
  H5Pset_deflate(dcpl1d, fCompressionLevel);

  fDsetInfo =
      H5Dcreate2(fFile, "/events/info", fEvtType, space1d, H5P_DEFAULT, dcpl1d, H5P_DEFAULT);
  fDsetIndex = H5Dcreate2(fFile, "/events/index", H5T_NATIVE_ULLONG, space1d, H5P_DEFAULT, dcpl1d,
                          H5P_DEFAULT);
  fDsetChs = H5Dcreate2(fFile, "/events/chs", fChType, space1d, H5P_DEFAULT, dcpl1d, H5P_DEFAULT);

  H5Pclose(dcpl1d);
  H5Sclose(space1d);

  // Create extendable 2D datasets for dual waveforms (phonon & photon)
  hsize_t dims2[2] = {0, static_cast<hsize_t>(fNDP)};
  hsize_t maxdims2[2] = {H5S_UNLIMITED, static_cast<hsize_t>(fNDP)};
  hid_t space2d = H5Screate_simple(2, dims2, maxdims2);

  hid_t dcpl2d = H5Pcreate(H5P_DATASET_CREATE);

  // CRITICAL TUNING: Chunk size along the event axis must be small to avoid huge memory allocation.
  // E.g., chunk[0]=2 means 2 * 500,000 * 2 bytes = ~2 MB per chunk block.
  hsize_t chunk2d[2] = {2, static_cast<hsize_t>(fNDP)};
  H5Pset_chunk(dcpl2d, 2, chunk2d);
  H5Pset_deflate(dcpl2d, fCompressionLevel);

  fDsetPhonon = H5Dcreate2(fFile, "/events/phonon", H5T_NATIVE_USHORT, space2d, H5P_DEFAULT, dcpl2d,
                           H5P_DEFAULT);
  fDsetPhoton = H5Dcreate2(fFile, "/events/photon", H5T_NATIVE_USHORT, space2d, H5P_DEFAULT, dcpl2d,
                           H5P_DEFAULT);

  H5Pclose(dcpl2d);
  H5Sclose(space2d);

  fMemSize = 0;
  fTotalEvents = 0;
  fTotalCrystals = 0;

  fEvtBuf.clear();
  fChBuf.clear();
  fPhononBuf.clear();
  fPhotonBuf.clear();
  fBufCount = 0;
  fBufBytesUsed = 0;
}

int H5AMOREEvent::GetNDP()
{
  if (fWriteTag || fNDP > 0) { return fNDP; }

  if (fFileSpacePhonon >= 0) {
    hsize_t dims[2] = {0, 0};
    H5Sget_simple_extent_dims(fFileSpacePhonon, dims, nullptr);
    fNDP = static_cast<int>(dims[1]);
    return fNDP;
  }

  hid_t fid = H5I_INVALID_HID;
  int dummy_evt = 0;

  if (fChain && fChain->GetNFile() > 0) { fid = fChain->GetFileId(0, dummy_evt); }
  else if (fFile >= 0) {
    fid = fFile;
  }

  if (fid < 0) { return 0; }

  hid_t dset = H5Dopen2(fid, "/events/phonon", H5P_DEFAULT);
  if (dset < 0) { return 0; }

  hid_t space = H5Dget_space(dset);
  hsize_t dims[2] = {0, 0};
  H5Sget_simple_extent_dims(space, dims, nullptr);
  H5Sclose(space);
  H5Dclose(dset);

  fNDP = (dims[1] > 0 && dims[1] <= static_cast<hsize_t>(kH5AMORENDPMAX))
             ? static_cast<int>(dims[1])
             : 0;
  return fNDP;
}

herr_t H5AMOREEvent::FlushBuffer()
{
  const std::size_t nEvtBuf = fEvtBuf.size();
  if (nEvtBuf == 0) {
    fBufCount = 0;
    fBufBytesUsed = 0;
    return 0;
  }

  const std::size_t nChBuf = fChBuf.size();
  herr_t status = 0;

  // 1. Append events info
  {
    hsize_t old_dim[1] = {fTotalEvents};
    hsize_t new_dim[1] = {fTotalEvents + nEvtBuf};
    H5Dset_extent(fDsetInfo, new_dim);

    hid_t file_space = H5Dget_space(fDsetInfo);
    hsize_t offset[1] = {old_dim[0]};
    hsize_t count[1] = {nEvtBuf};

    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    hid_t mem_space = H5Screate_simple(1, count, nullptr);
    status = H5Dwrite(fDsetInfo, fEvtType, mem_space, file_space, H5P_DEFAULT, fEvtBuf.data());

    H5Sclose(mem_space);
    H5Sclose(file_space);
    if (status < 0) return status;
  }

  // 2. Append event index
  {
    std::vector<std::uint64_t> indexBuf(nEvtBuf);
    std::uint64_t base = fTotalCrystals;
    std::uint64_t running = 0;
    for (std::size_t i = 0; i < nEvtBuf; ++i) {
      indexBuf[i] = base + running;
      running += static_cast<std::uint64_t>(fEvtBuf[i].nhit);
    }

    hsize_t old_dim[1] = {fTotalEvents};
    hsize_t new_dim[1] = {fTotalEvents + nEvtBuf};
    H5Dset_extent(fDsetIndex, new_dim);

    hid_t file_space = H5Dget_space(fDsetIndex);
    hsize_t offset[1] = {old_dim[0]};
    hsize_t count[1] = {nEvtBuf};

    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    hid_t mem_space = H5Screate_simple(1, count, nullptr);
    status = H5Dwrite(fDsetIndex, H5T_NATIVE_ULLONG, mem_space, file_space, H5P_DEFAULT,
                      indexBuf.data());

    H5Sclose(mem_space);
    H5Sclose(file_space);
    if (status < 0) return status;
  }

  // 3. Append channels, phonon, and photon waves
  if (nChBuf > 0) {
    hsize_t old_dim[1] = {fTotalCrystals};
    hsize_t new_dim[1] = {fTotalCrystals + nChBuf};
    H5Dset_extent(fDsetChs, new_dim);

    hid_t file_space = H5Dget_space(fDsetChs);
    hsize_t offset[1] = {old_dim[0]};
    hsize_t count[1] = {nChBuf};

    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    hid_t mem_space = H5Screate_simple(1, count, nullptr);
    status = H5Dwrite(fDsetChs, fChType, mem_space, file_space, H5P_DEFAULT, fChBuf.data());

    H5Sclose(mem_space);
    H5Sclose(file_space);
    if (status < 0) return status;

    // Expand 2D datasets
    hsize_t old_dims[2] = {fTotalCrystals, static_cast<hsize_t>(fNDP)};
    hsize_t new_dims[2] = {fTotalCrystals + nChBuf, static_cast<hsize_t>(fNDP)};
    H5Dset_extent(fDsetPhonon, new_dims);
    H5Dset_extent(fDsetPhoton, new_dims);

    hsize_t offset_wave[2] = {old_dims[0], 0};
    hsize_t count_wave[2] = {nChBuf, static_cast<hsize_t>(fNDP)};
    hid_t mem_space_wave = H5Screate_simple(2, count_wave, nullptr);

    // Write Phonon
    hid_t file_space_pn = H5Dget_space(fDsetPhonon);
    H5Sselect_hyperslab(file_space_pn, H5S_SELECT_SET, offset_wave, nullptr, count_wave, nullptr);
    status = H5Dwrite(fDsetPhonon, H5T_NATIVE_USHORT, mem_space_wave, file_space_pn, H5P_DEFAULT,
                      fPhononBuf.data());
    H5Sclose(file_space_pn);
    if (status < 0) {
      H5Sclose(mem_space_wave);
      return status;
    }

    // Write Photon
    hid_t file_space_pt = H5Dget_space(fDsetPhoton);
    H5Sselect_hyperslab(file_space_pt, H5S_SELECT_SET, offset_wave, nullptr, count_wave, nullptr);
    status = H5Dwrite(fDsetPhoton, H5T_NATIVE_USHORT, mem_space_wave, file_space_pt, H5P_DEFAULT,
                      fPhotonBuf.data());
    H5Sclose(file_space_pt);

    H5Sclose(mem_space_wave);
    if (status < 0) return status;
  }

  fTotalEvents += static_cast<std::uint64_t>(nEvtBuf);
  fTotalCrystals += static_cast<std::uint64_t>(nChBuf);

  fEvtBuf.clear();
  fChBuf.clear();
  fPhononBuf.clear();
  fPhotonBuf.clear();
  fBufCount = 0;
  fBufBytesUsed = 0;

  return status;
}

void H5AMOREEvent::Close()
{
  if (fWriteTag) { FlushBuffer(); }

  if (fFileSpaceInfo >= 0) {
    H5Sclose(fFileSpaceInfo);
    fFileSpaceInfo = H5I_INVALID_HID;
  }
  if (fFileSpaceIndex >= 0) {
    H5Sclose(fFileSpaceIndex);
    fFileSpaceIndex = H5I_INVALID_HID;
  }
  if (fFileSpaceChs >= 0) {
    H5Sclose(fFileSpaceChs);
    fFileSpaceChs = H5I_INVALID_HID;
  }
  if (fFileSpacePhonon >= 0) {
    H5Sclose(fFileSpacePhonon);
    fFileSpacePhonon = H5I_INVALID_HID;
  }
  if (fFileSpacePhoton >= 0) {
    H5Sclose(fFileSpacePhoton);
    fFileSpacePhoton = H5I_INVALID_HID;
  }

  if (fDsetInfo >= 0) {
    H5Dclose(fDsetInfo);
    fDsetInfo = H5I_INVALID_HID;
  }
  if (fDsetIndex >= 0) {
    H5Dclose(fDsetIndex);
    fDsetIndex = H5I_INVALID_HID;
  }
  if (fDsetChs >= 0) {
    H5Dclose(fDsetChs);
    fDsetChs = H5I_INVALID_HID;
  }
  if (fDsetPhonon >= 0) {
    H5Dclose(fDsetPhonon);
    fDsetPhonon = H5I_INVALID_HID;
  }
  if (fDsetPhoton >= 0) {
    H5Dclose(fDsetPhoton);
    fDsetPhoton = H5I_INVALID_HID;
  }

  fCurrentReadFid = H5I_INVALID_HID;

  fReadBufInfo.clear();
  fReadBufIndex.clear();
  fPrefetchChs.clear();
  fPrefetchPhonon.clear();
  fPrefetchPhoton.clear();

  if (fEvtType >= 0) {
    H5Tclose(fEvtType);
    fEvtType = H5I_INVALID_HID;
  }
  if (fChType >= 0) {
    H5Tclose(fChType);
    fChType = H5I_INVALID_HID;
  }

  CloseSubRun();
}

herr_t H5AMOREEvent::AppendEvent(const EventInfo_t & info, const std::vector<Crystal_t> & data)
{
  if (!fWriteTag) { return -1; }
  if (fNDP <= 0 || fNDP > kH5AMORENDPMAX) {
    Error("AppendEvent", "Invalid NDP: %d", fNDP);
    return -1;
  }

  const std::uint16_t nhit = static_cast<std::uint16_t>(data.size());

  EventInfo_t info_local = info;
  info_local.nhit = nhit;
  fEvtBuf.push_back(info_local);

  std::size_t addBytes = sizeof(EventInfo_t);
  // Add metadata size + dual waveform payload size per hit
  addBytes +=
      static_cast<std::size_t>(nhit) *
      (sizeof(CrystalHeader_t) + 2 * static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t));

  for (std::size_t i = 0; i < nhit; ++i) {
    CrystalHeader_t h{};
    h.id = data[i].id;
    h.ttime = data[i].ttime;
    fChBuf.push_back(h);

    const std::uint16_t * src_pn = data[i].phonon;
    fPhononBuf.insert(fPhononBuf.end(), src_pn, src_pn + static_cast<std::size_t>(fNDP));

    const std::uint16_t * src_pt = data[i].photon;
    fPhotonBuf.insert(fPhotonBuf.end(), src_pt, src_pt + static_cast<std::size_t>(fNDP));
  }

  UpdateSubRun(info.tnum);

  fMemSize += static_cast<hsize_t>(addBytes);
  fBufCount += 1;
  fBufBytesUsed += addBytes;

  if ((fBufCap > 0 && fBufCount >= fBufCap) ||
      (fBufMaxBytes > 0 && fBufBytesUsed >= fBufMaxBytes)) {
    return FlushBuffer();
  }

  return 0;
}

herr_t H5AMOREEvent::ReadEvent(int n)
{
  int evtno = n;
  bool file_changed = false;
  hid_t fid = H5I_INVALID_HID;

  if (fChain && fChain->GetNFile() > 0) { fid = fChain->GetFileId(n, evtno, &file_changed); }
  else {
    fid = fFile;
    if (fCurrentReadFid != fid) { file_changed = true; }
  }

  if (fid < 0) return -1;

  if (file_changed) {
    if (fFileSpaceInfo >= 0) { H5Sclose(fFileSpaceInfo); }
    if (fFileSpaceIndex >= 0) { H5Sclose(fFileSpaceIndex); }
    if (fFileSpaceChs >= 0) { H5Sclose(fFileSpaceChs); }
    if (fFileSpacePhonon >= 0) { H5Sclose(fFileSpacePhonon); }
    if (fFileSpacePhoton >= 0) { H5Sclose(fFileSpacePhoton); }

    if (fDsetInfo >= 0) { H5Dclose(fDsetInfo); }
    if (fDsetIndex >= 0) { H5Dclose(fDsetIndex); }
    if (fDsetChs >= 0) { H5Dclose(fDsetChs); }
    if (fDsetPhonon >= 0) { H5Dclose(fDsetPhonon); }
    if (fDsetPhoton >= 0) { H5Dclose(fDsetPhoton); }

    hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk_cache(dapl, 10007, 256 * 1024 * 1024,
                       1.0); // Increase cache for large AMORE arrays

    fDsetInfo = H5Dopen2(fid, "/events/info", H5P_DEFAULT);
    fDsetIndex = H5Dopen2(fid, "/events/index", H5P_DEFAULT);
    fDsetChs = H5Dopen2(fid, "/events/chs", H5P_DEFAULT);
    fDsetPhonon = H5Dopen2(fid, "/events/phonon", dapl);
    fDsetPhoton = H5Dopen2(fid, "/events/photon", dapl);
    H5Pclose(dapl);

    if (fDsetInfo < 0 || fDsetIndex < 0 || fDsetChs < 0 || fDsetPhonon < 0 || fDsetPhoton < 0)
      return -1;

    fFileSpaceInfo = H5Dget_space(fDsetInfo);
    fFileSpaceIndex = H5Dget_space(fDsetIndex);
    fFileSpaceChs = H5Dget_space(fDsetChs);
    fFileSpacePhonon = H5Dget_space(fDsetPhonon);
    fFileSpacePhoton = H5Dget_space(fDsetPhoton);

    fCurrentReadFid = fid;
    fNDP = GetNDP();
    fReadBufStart = -1;
    fReadBufSize = 0;
  }

  if (evtno < fReadBufStart || evtno >= fReadBufStart + fReadBufSize) {
    hsize_t dims[1];
    H5Sget_simple_extent_dims(fFileSpaceInfo, dims, nullptr);
    int total_events = static_cast<int>(dims[0]);

    std::size_t bytes_per_event = sizeof(EventInfo_t) + sizeof(std::uint64_t) +
                                  4 * (sizeof(CrystalHeader_t) + 2 * fNDP * sizeof(std::uint16_t));
    int max_safe_events = static_cast<int>(fBufMaxBytes / bytes_per_event);
    if (max_safe_events <= 0) max_safe_events = 1; // At least prefetch 1 to prevent deadlock

    int fetch_size = std::min(max_safe_events, total_events - evtno);

    if (fetch_size > 0) {
      fReadBufInfo.resize(fetch_size);
      fReadBufIndex.resize(fetch_size);

      hsize_t offset_blk[1] = {static_cast<hsize_t>(evtno)};
      hsize_t count_blk[1] = {static_cast<hsize_t>(fetch_size)};
      hid_t mem_space_blk = H5Screate_simple(1, count_blk, nullptr);

      H5Sselect_hyperslab(fFileSpaceInfo, H5S_SELECT_SET, offset_blk, nullptr, count_blk, nullptr);
      H5Dread(fDsetInfo, fEvtType, mem_space_blk, fFileSpaceInfo, H5P_DEFAULT, fReadBufInfo.data());

      H5Sselect_hyperslab(fFileSpaceIndex, H5S_SELECT_SET, offset_blk, nullptr, count_blk, nullptr);
      H5Dread(fDsetIndex, H5T_NATIVE_ULLONG, mem_space_blk, fFileSpaceIndex, H5P_DEFAULT,
              fReadBufIndex.data());
      H5Sclose(mem_space_blk);

      std::uint64_t first_ch = fReadBufIndex[0];
      std::uint64_t total_chs =
          fReadBufIndex[fetch_size - 1] + fReadBufInfo[fetch_size - 1].nhit - first_ch;

      if (total_chs > 0) {
        fPrefetchChs.resize(total_chs);
        fPrefetchPhonon.resize(total_chs * fNDP);
        fPrefetchPhoton.resize(total_chs * fNDP);

        hsize_t offset_ch[1] = {first_ch};
        hsize_t count_ch[1] = {total_chs};
        H5Sselect_hyperslab(fFileSpaceChs, H5S_SELECT_SET, offset_ch, nullptr, count_ch, nullptr);
        hid_t mem_space_chs = H5Screate_simple(1, count_ch, nullptr);
        H5Dread(fDsetChs, fChType, mem_space_chs, fFileSpaceChs, H5P_DEFAULT, fPrefetchChs.data());
        H5Sclose(mem_space_chs);

        hsize_t offset_wave[2] = {first_ch, 0};
        hsize_t count_wave[2] = {total_chs, static_cast<hsize_t>(fNDP)};
        hid_t mem_space_wave = H5Screate_simple(2, count_wave, nullptr);

        H5Sselect_hyperslab(fFileSpacePhonon, H5S_SELECT_SET, offset_wave, nullptr, count_wave,
                            nullptr);
        H5Dread(fDsetPhonon, H5T_NATIVE_USHORT, mem_space_wave, fFileSpacePhonon, H5P_DEFAULT,
                fPrefetchPhonon.data());

        H5Sselect_hyperslab(fFileSpacePhoton, H5S_SELECT_SET, offset_wave, nullptr, count_wave,
                            nullptr);
        H5Dread(fDsetPhoton, H5T_NATIVE_USHORT, mem_space_wave, fFileSpacePhoton, H5P_DEFAULT,
                fPrefetchPhoton.data());

        H5Sclose(mem_space_wave);
      }

      fReadBufStart = evtno;
      fReadBufSize = fetch_size;
      fPrefetchChStart = first_ch;
    }
    else {
      return -1;
    }
  }

  int local_idx = evtno - fReadBufStart;
  fEvtInfo = fReadBufInfo[local_idx];
  std::uint64_t ch_offset = fReadBufIndex[local_idx];

  const std::uint16_t nhit = fEvtInfo.nhit;
  fDataBuf.resize(nhit);

  if (nhit > 0) {
    std::uint64_t local_ch_idx = ch_offset - fPrefetchChStart;

    for (std::size_t ich = 0; ich < nhit; ++ich) {
      fDataBuf[ich].id = fPrefetchChs[local_ch_idx + ich].id;
      fDataBuf[ich].ttime = fPrefetchChs[local_ch_idx + ich].ttime;

      std::uint16_t * dst_pn = fDataBuf[ich].phonon;
      const std::uint16_t * src_pn = &fPrefetchPhonon[(local_ch_idx + ich) * fNDP];
      std::memcpy(dst_pn, src_pn, static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t));

      std::uint16_t * dst_pt = fDataBuf[ich].photon;
      const std::uint16_t * src_pt = &fPrefetchPhoton[(local_ch_idx + ich) * fNDP];
      std::memcpy(dst_pt, src_pt, static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t));
    }
  }

  return 0;
}
#include <algorithm>
#include <cstring>

#include "AMOREHDF5/H5AMOREHit.hh"

ClassImp(H5AMOREHit)

H5AMOREHit::H5AMOREHit()
  : AbsH5Hit()
{
}

H5AMOREHit::~H5AMOREHit() {}

void H5AMOREHit::Open()
{
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

  // Create /hits group for Self Trigger stream
  {
    htri_t gexists = H5Lexists(fFile, "/hits", H5P_DEFAULT);
    if (gexists < 0) {
      Error("Open", "H5Lexists(/hits) failed");
      return;
    }
    if (gexists == 0) {
      hid_t grp_hits = H5Gcreate2(fFile, "/hits", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (grp_hits < 0) {
        Error("Open", "Failed to create group /hits");
        return;
      }
      H5Gclose(grp_hits);
    }
  }

  // Create extendable dataset: /hits/chs (flattened hit metadata)
  {
    hsize_t dims[1] = {0};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hid_t space = H5Screate_simple(1, dims, maxdims);

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk = 1024;
    H5Pset_chunk(dcpl, 1, &chunk);
    H5Pset_deflate(dcpl, fCompressionLevel);

    fDsetChs = H5Dcreate2(fFile, "/hits/chs", fChType, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    H5Pclose(dcpl);
    H5Sclose(space);

    if (fDsetChs < 0) {
      Error("Open", "Failed to create dataset /hits/chs");
      return;
    }
  }

  // Create extendable datasets: /hits/phonon & /hits/photon
  {
    hsize_t dims[2] = {0, static_cast<hsize_t>(fNDP)};
    hsize_t maxdims[2] = {H5S_UNLIMITED, static_cast<hsize_t>(fNDP)};
    hid_t space = H5Screate_simple(2, dims, maxdims);

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);

    // CRITICAL: Extremely small chunk size along the hit dimension to manage massive memory
    // footprint
    hsize_t chunk[2] = {2, static_cast<hsize_t>(fNDP)};
    H5Pset_chunk(dcpl, 2, chunk);
    H5Pset_deflate(dcpl, fCompressionLevel);

    fDsetPhonon =
        H5Dcreate2(fFile, "/hits/phonon", H5T_NATIVE_USHORT, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    fDsetPhoton =
        H5Dcreate2(fFile, "/hits/photon", H5T_NATIVE_USHORT, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    H5Pclose(dcpl);
    H5Sclose(space);

    if (fDsetPhonon < 0 || fDsetPhoton < 0) {
      Error("Open", "Failed to create dual waveform datasets");
      return;
    }
  }

  fMemSize = 0;
  fTotalHits = 0;

  fChBuf.clear();
  fPhononBuf.clear();
  fPhotonBuf.clear();
  fBufCount = 0;
  fBufBytesUsed = 0;
}

int H5AMOREHit::GetNDP()
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

  if (fid < 0) return 0;

  hid_t dset = H5Dopen2(fid, "/hits/phonon", H5P_DEFAULT);
  if (dset < 0) return 0;

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

herr_t H5AMOREHit::FlushBuffer()
{
  const std::size_t nHitBuf = fChBuf.size();
  if (nHitBuf == 0) {
    fBufCount = 0;
    fBufBytesUsed = 0;
    return 0;
  }

  herr_t status = 0;

  // Append hit headers to /hits/chs
  {
    hsize_t old_dim[1] = {fTotalHits};
    hsize_t new_dim[1] = {fTotalHits + nHitBuf};
    H5Dset_extent(fDsetChs, new_dim);

    hid_t file_space = H5Dget_space(fDsetChs);
    hsize_t offset[1] = {old_dim[0]};
    hsize_t count[1] = {nHitBuf};
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    hid_t mem_space = H5Screate_simple(1, count, nullptr);
    status = H5Dwrite(fDsetChs, fChType, mem_space, file_space, H5P_DEFAULT, fChBuf.data());

    H5Sclose(mem_space);
    H5Sclose(file_space);
    if (status < 0) return status;
  }

  // Append dual waveforms to /hits/phonon & /hits/photon
  {
    hsize_t old_dims[2] = {fTotalHits, static_cast<hsize_t>(fNDP)};
    hsize_t new_dims[2] = {fTotalHits + nHitBuf, static_cast<hsize_t>(fNDP)};
    H5Dset_extent(fDsetPhonon, new_dims);
    H5Dset_extent(fDsetPhoton, new_dims);

    hsize_t offset[2] = {old_dims[0], 0};
    hsize_t count[2] = {nHitBuf, static_cast<hsize_t>(fNDP)};
    hid_t mem_space = H5Screate_simple(2, count, nullptr);

    // Write Phonon
    hid_t file_space_pn = H5Dget_space(fDsetPhonon);
    H5Sselect_hyperslab(file_space_pn, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite(fDsetPhonon, H5T_NATIVE_USHORT, mem_space, file_space_pn, H5P_DEFAULT,
                      fPhononBuf.data());
    H5Sclose(file_space_pn);
    if (status < 0) {
      H5Sclose(mem_space);
      return status;
    }

    // Write Photon
    hid_t file_space_pt = H5Dget_space(fDsetPhoton);
    H5Sselect_hyperslab(file_space_pt, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite(fDsetPhoton, H5T_NATIVE_USHORT, mem_space, file_space_pt, H5P_DEFAULT,
                      fPhotonBuf.data());
    H5Sclose(file_space_pt);

    H5Sclose(mem_space);
    if (status < 0) return status;
  }

  fTotalHits += static_cast<std::uint64_t>(nHitBuf);

  fChBuf.clear();
  fPhononBuf.clear();
  fPhotonBuf.clear();
  fBufCount = 0;
  fBufBytesUsed = 0;

  return status;
}

void H5AMOREHit::Close()
{
  if (fWriteTag) { FlushBuffer(); }

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

  fPrefetchChs.clear();
  fPrefetchPhonon.clear();
  fPrefetchPhoton.clear();

  if (fChType >= 0) {
    H5Tclose(fChType);
    fChType = H5I_INVALID_HID;
  }

  CloseSubRun();
}

herr_t H5AMOREHit::AppendHit(const Crystal_t & hit)
{
  if (!fWriteTag) return -1;
  if (fNDP <= 0 || fNDP > kH5AMORENDPMAX) return -1;

  CrystalHeader_t header;
  header.id = hit.id;
  header.ttime = hit.ttime;

  fChBuf.push_back(header);
  fPhononBuf.insert(fPhononBuf.end(), hit.phonon, hit.phonon + static_cast<std::size_t>(fNDP));
  fPhotonBuf.insert(fPhotonBuf.end(), hit.photon, hit.photon + static_cast<std::size_t>(fNDP));

  std::size_t addBytes =
      sizeof(CrystalHeader_t) + 2 * static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t);
  fMemSize += static_cast<hsize_t>(addBytes);
  fBufCount += 1;
  fBufBytesUsed += addBytes;

  // Update SubRun Tracker
  UpdateSubRun(header.ttime);

  if ((fBufCap > 0 && fBufCount >= fBufCap) ||
      (fBufMaxBytes > 0 && fBufBytesUsed >= fBufMaxBytes)) {
    return FlushBuffer();
  }

  return 0;
}

herr_t H5AMOREHit::ReadHit(std::uint64_t n)
{
  std::uint64_t local_hit_no = n;
  hid_t fid = H5I_INVALID_HID;

  if (fChain && fChain->GetNFile() > 0) {
    int dummy = 0;
    fid = fChain->GetFileId(0, dummy);
  }
  else {
    fid = fFile;
  }

  if (fid < 0) return -1;

  if (fCurrentReadFid != fid) {
    if (fFileSpaceChs >= 0) H5Sclose(fFileSpaceChs);
    if (fFileSpacePhonon >= 0) H5Sclose(fFileSpacePhonon);
    if (fFileSpacePhoton >= 0) H5Sclose(fFileSpacePhoton);
    if (fDsetChs >= 0) H5Dclose(fDsetChs);
    if (fDsetPhonon >= 0) H5Dclose(fDsetPhonon);
    if (fDsetPhoton >= 0) H5Dclose(fDsetPhoton);

    hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk_cache(dapl, 10007, 256 * 1024 * 1024, 1.0);

    fDsetChs = H5Dopen2(fid, "/hits/chs", H5P_DEFAULT);
    fDsetPhonon = H5Dopen2(fid, "/hits/phonon", dapl);
    fDsetPhoton = H5Dopen2(fid, "/hits/photon", dapl);
    H5Pclose(dapl);

    if (fDsetChs < 0 || fDsetPhonon < 0 || fDsetPhoton < 0) return -1;

    fFileSpaceChs = H5Dget_space(fDsetChs);
    fFileSpacePhonon = H5Dget_space(fDsetPhonon);
    fFileSpacePhoton = H5Dget_space(fDsetPhoton);

    fCurrentReadFid = fid;
    fNDP = GetNDP();
    fReadBufStart = 0;
    fReadBufSize = 0;
  }

  if (local_hit_no < fReadBufStart || local_hit_no >= fReadBufStart + fReadBufSize) {
    hsize_t dims[1];
    H5Sget_simple_extent_dims(fFileSpaceChs, dims, nullptr);
    std::uint64_t total_hits_in_file = dims[0];

    if (local_hit_no >= total_hits_in_file) return -1;

    std::size_t bytes_per_hit =
        sizeof(CrystalHeader_t) + 2 * static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t);
    if (bytes_per_hit == 0) return -1;

    std::uint64_t max_safe_hits = fBufMaxBytes / bytes_per_hit;
    if (max_safe_hits == 0) max_safe_hits = 1;

    std::uint64_t remaining = total_hits_in_file - local_hit_no;
    std::uint64_t fetch_size = std::min(max_safe_hits, remaining);

    if (fetch_size > 0) {
      fPrefetchChs.resize(fetch_size);
      fPrefetchPhonon.resize(fetch_size * fNDP);
      fPrefetchPhoton.resize(fetch_size * fNDP);

      hsize_t offset_ch[1] = {local_hit_no};
      hsize_t count_ch[1] = {fetch_size};
      hid_t mem_space_chs = H5Screate_simple(1, count_ch, nullptr);
      H5Sselect_hyperslab(fFileSpaceChs, H5S_SELECT_SET, offset_ch, nullptr, count_ch, nullptr);
      H5Dread(fDsetChs, fChType, mem_space_chs, fFileSpaceChs, H5P_DEFAULT, fPrefetchChs.data());
      H5Sclose(mem_space_chs);

      hsize_t offset_wave[2] = {local_hit_no, 0};
      hsize_t count_wave[2] = {fetch_size, static_cast<hsize_t>(fNDP)};
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

      fReadBufStart = local_hit_no;
      fReadBufSize = fetch_size;
    }
  }

  std::uint64_t local_idx = local_hit_no - fReadBufStart;

  fCurrentHit.id = fPrefetchChs[local_idx].id;
  fCurrentHit.ttime = fPrefetchChs[local_idx].ttime;

  std::memcpy(fCurrentHit.phonon, &fPrefetchPhonon[local_idx * fNDP],
              static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t));
  std::memcpy(fCurrentHit.photon, &fPrefetchPhoton[local_idx * fNDP],
              static_cast<std::size_t>(fNDP) * sizeof(std::uint16_t));

  return 0;
}
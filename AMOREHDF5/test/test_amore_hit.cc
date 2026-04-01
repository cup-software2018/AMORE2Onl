#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "AMOREHDF5/H5AMOREHit.hh"
#include "HDF5Utils/H5DataReader.hh"
#include "HDF5Utils/H5DataWriter.hh"

// -----------------------------------------------------------------------------
// Write function for AMORE Hits (Self Trigger)
// -----------------------------------------------------------------------------
void WriteAMOREHit(const char * filename, int n_hits, int ndp)
{
  // Initialize number generator for variable crystal IDs and ADC values
  std::mt19937 rng(54321);                               // Fixed seed for reproducibility
  std::uniform_int_distribution<int> dist_id(0, 39);     // Assume 40 crystals in total
  std::uniform_int_distribution<int> dist_adc(0, 65535); // 16-bit range

  std::cout << "========== [WRITE] Start writing AMORE Hits ==========\n";

  // Setup DataWriter
  H5DataWriter writer(filename);
  writer.SetCompressionLevel(4);

  H5AMOREHit * wHit = new H5AMOREHit();
  wHit->SetNDP(ndp);
  wHit->SetBufferCapacity(100); // Flush to file safely

  writer.SetData(wHit);
  writer.SetSubrun(999);

  // Start the timer
  auto start_time = std::chrono::high_resolution_clock::now();

  writer.Open();

  // Generate and append independent hit data
  for (int i = 0; i < n_hits; ++i) {
    Crystal_t current_hit;
    current_hit.id = static_cast<std::uint16_t>(dist_id(rng));
    current_hit.ttime = static_cast<std::uint64_t>(i) * 5000ULL; // Simulated trigger time

    // Fill dual waveform data (Phonon & Photon)
    for (int s = 0; s < ndp; ++s) {
      current_hit.phonon[s] = static_cast<std::uint16_t>(dist_adc(rng));
      current_hit.photon[s] = static_cast<std::uint16_t>(dist_adc(rng));
    }

    wHit->AppendHit(current_hit);
  }

  // Get file size before closing the file
  double file_size_mb = static_cast<double>(writer.GetFileSize()) / (1024.0 * 1024.0);

  writer.PrintStats();
  writer.Close();

  // Stop the timer
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  double elapsed_sec = elapsed.count();

  // ---------------------------------------------------------------------------
  // Calculate Statistics (Speed and Compression)
  // ---------------------------------------------------------------------------

  // Calculate raw uncompressed bytes
  double raw_header_bytes = static_cast<double>(n_hits * sizeof(CrystalHeader_t));
  // Multiply by 2 because AMORE has BOTH phonon and photon arrays
  double raw_wave_bytes = static_cast<double>(n_hits * ndp * 2 * sizeof(std::uint16_t));

  double uncompressed_mb = (raw_header_bytes + raw_wave_bytes) / (1024.0 * 1024.0);

  // Avoid division by zero if file creation failed
  double compression_ratio = (file_size_mb > 0.0) ? (uncompressed_mb / file_size_mb) : 0.0;

  double speed_mb_per_sec = file_size_mb / elapsed_sec;
  double speed_hits_per_sec = n_hits / elapsed_sec;

  // Print detailed write statistics
  std::cout << "\n--- Write Statistics Summary ---\n";
  std::cout << "Total Hits Written       : " << n_hits << "\n";
  std::cout << "Data Points Per Hit      : " << ndp << " (x2 arrays)\n";
  std::cout << "--------------------------------\n";
  std::cout << "Uncompressed Data Size   : " << std::fixed << std::setprecision(2)
            << uncompressed_mb << " MB\n";
  std::cout << "Compressed File Size     : " << std::fixed << std::setprecision(2) << file_size_mb
            << " MB\n";
  std::cout << "Compression Ratio        : " << std::fixed << std::setprecision(2)
            << compression_ratio << "x\n";
  std::cout << "--------------------------------\n";
  std::cout << "Elapsed Time             : " << std::fixed << std::setprecision(4) << elapsed_sec
            << " seconds\n";
  std::cout << "Write Speed (Disk IO)    : " << std::fixed << std::setprecision(2)
            << speed_mb_per_sec << " MB/s\n";
  std::cout << "Write Speed (Hits)       : " << std::fixed << std::setprecision(2)
            << speed_hits_per_sec << " Hits/s\n";
  std::cout << "--------------------------------\n\n";
}

// -----------------------------------------------------------------------------
// Read function for AMORE Hits
// -----------------------------------------------------------------------------
void ReadAMOREHit(const char * filename)
{
  std::cout << "========== [READ] Read and verify AMORE Hits ==========\n";

  // Setup DataReader
  H5DataReader reader;
  H5AMOREHit * rHit = new H5AMOREHit();

  reader.Add(filename);
  reader.SetData(rHit);

  if (!reader.Open()) {
    std::cerr << "Failed to open the file.\n";
    return;
  }

  int total_entries = reader.GetEntries();
  std::cout << "Total hits to read: " << total_entries << "\n\n";

  // Read and verify data
  int check_interval = std::max(1, total_entries / 10); // Print roughly 10 checkpoints

  for (int i = 0; i < total_entries; ++i) {
    if (rHit->ReadHit(i) < 0) {
      std::cerr << "Failed to read hit " << i << "!\n";
      break;
    }

    if (i % check_interval == 0 || i == total_entries - 1) {
      const Crystal_t & hit_data = rHit->GetHit();

      std::cout << "[Hit " << i << "] "
                << "Crystal ID: " << hit_data.id << " | TTime: " << hit_data.ttime << "\n"
                << "  -> Phonon[0]: " << hit_data.phonon[0]
                << " | Photon[0]: " << hit_data.photon[0] << "\n";
      std::cout << "--------------------------------------------------\n";
    }
  }

  reader.Close();
  std::cout << "AMORE Hit read test completed successfully.\n";
}

// -----------------------------------------------------------------------------
// Main execution
// -----------------------------------------------------------------------------
int main()
{
  const char * filename = "test_amore_hits.h5";

  // Adjusted parameters for AMORE's massive data size
  // 2500 hits * 20,000 NDP * 2 (arrays) * 2 bytes = ~200 MB uncompressed
  const int n_hits = 2500;
  const int ndp = 20000;

  WriteAMOREHit(filename, n_hits, ndp);
  ReadAMOREHit(filename);

  return 0;
}
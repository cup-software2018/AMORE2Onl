#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "AMOREHDF5/H5AMOREEvent.hh"
#include "HDF5Utils/H5DataReader.hh"
#include "HDF5Utils/H5DataWriter.hh"

// -----------------------------------------------------------------------------
// Write function for AMORE Events
// -----------------------------------------------------------------------------
void WriteAMOREEvent(const char * filename, int n_events, int max_crystals, int ndp)
{
  // Initialize number generator for variable crystals and ADC values
  std::mt19937 rng(12345); // Fixed seed for reproducibility
  std::uniform_int_distribution<int> dist_nhit(1, max_crystals);
  std::uniform_int_distribution<int> dist_adc(0, 65535); // 16-bit range for AMORE

  std::cout << "========== [WRITE] Start writing AMORE events ==========\n";

  // Setup DataWriter
  H5DataWriter writer(filename);
  writer.SetCompressionLevel(4);

  H5AMOREEvent * wEvent = new H5AMOREEvent();
  wEvent->SetNDP(ndp);
  wEvent->SetBufferCapacity(100); // Flush to file safely (memory managed by fBufMaxBytes as well)

  writer.SetData(wEvent);
  writer.SetSubrun(999);

  // Start the timer
  auto start_time = std::chrono::high_resolution_clock::now();

  writer.Open();

  // Generate and append event data
  int total_crystals_written = 0;

  for (int i = 0; i < n_events; ++i) {
    EventInfo_t info;
    info.tnum = i;
    info.ttime = static_cast<std::uint64_t>(i) * 10000ULL;

    // Determine the number of triggered crystals for the current event
    int current_nhit = dist_nhit(rng);
    total_crystals_written += current_nhit;

    std::vector<Crystal_t> crystals(current_nhit);
    for (int ch = 0; ch < current_nhit; ++ch) {
      crystals[ch].id = ch;
      crystals[ch].ttime = info.ttime + ch;

      // Fill dual waveform data (Phonon & Photon)
      for (int s = 0; s < ndp; ++s) {
        crystals[ch].phonon[s] = static_cast<std::uint16_t>(dist_adc(rng));
        crystals[ch].photon[s] = static_cast<std::uint16_t>(dist_adc(rng));
      }
    }

    wEvent->AppendEvent(info, crystals);
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

  // Calculate raw uncompressed bytes based on AMORE EDM struct sizes
  double raw_event_bytes = static_cast<double>(n_events * sizeof(EventInfo_t));
  double raw_header_bytes = static_cast<double>(total_crystals_written * sizeof(CrystalHeader_t));

  // Multiply by 2 because AMORE has BOTH phonon and photon arrays
  double raw_wave_bytes =
      static_cast<double>(total_crystals_written * ndp * 2 * sizeof(std::uint16_t));

  double uncompressed_mb =
      (raw_event_bytes + raw_header_bytes + raw_wave_bytes) / (1024.0 * 1024.0);

  // Avoid division by zero if file creation failed
  double compression_ratio = (file_size_mb > 0.0) ? (uncompressed_mb / file_size_mb) : 0.0;

  double speed_mb_per_sec = file_size_mb / elapsed_sec;
  double speed_events_per_sec = n_events / elapsed_sec;

  // Print detailed write statistics
  std::cout << "\n--- Write Statistics Summary ---\n";
  std::cout << "Total Events Written     : " << n_events << "\n";
  std::cout << "Total Crystals Written   : " << total_crystals_written << "\n";
  std::cout << "Data Points Per Crystal  : " << ndp << " (x2 arrays)\n";
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
  std::cout << "Write Speed (Events)     : " << std::fixed << std::setprecision(2)
            << speed_events_per_sec << " Events/s\n";
  std::cout << "--------------------------------\n\n";
}

// -----------------------------------------------------------------------------
// Read function for AMORE Events
// -----------------------------------------------------------------------------
void ReadAMOREEvent(const char * filename)
{
  std::cout << "========== [READ] Read and verify AMORE events ==========\n";

  // Setup DataReader
  H5DataReader reader;
  H5AMOREEvent * rEvent = new H5AMOREEvent();

  reader.Add(filename);
  reader.SetData(rEvent);

  if (!reader.Open()) {
    std::cerr << "Failed to open the file.\n";
    return;
  }

  int total_entries = reader.GetEntries();
  std::cout << "Total events to read: " << total_entries << "\n\n";

  // Read and verify data
  int check_interval = std::max(1, total_entries / 10); // Print roughly 10 checkpoints

  for (int i = 0; i < total_entries; ++i) {
    if (rEvent->ReadEvent(i) < 0) {
      std::cerr << "Failed to read event " << i << "!\n";
      break;
    }

    if (i % check_interval == 0 || i == total_entries - 1) {
      EventInfo_t info = rEvent->GetEventInfo();
      Crystal_t * data = rEvent->GetData();

      std::cout << "[Event " << i << "] TNum: " << info.tnum
                << ", Triggered Crystals (NHits): " << info.nhit << "\n";

      if (info.nhit > 0) {
        std::cout << "  -> Crystal " << data[0].id << " | Phonon[0]: " << data[0].phonon[0]
                  << " | Photon[0]: " << data[0].photon[0] << "\n";

        if (info.nhit > 1) {
          int last_idx = info.nhit - 1;
          std::cout << "  -> Crystal " << data[last_idx].id
                    << " | Phonon[0]: " << data[last_idx].phonon[0]
                    << " | Photon[0]: " << data[last_idx].photon[0] << "\n";
        }
      }
      std::cout << "--------------------------------------------------\n";
    }
  }

  reader.Close();
  std::cout << "AMORE Event read test completed successfully.\n";
}

// -----------------------------------------------------------------------------
// Main execution
// -----------------------------------------------------------------------------
int main()
{
  const char * filename = "test_amore_events.h5";

  // Adjusted parameters for AMORE's massive data size
  // 500 events * ~5 crystals * 20,000 NDP * 2 (arrays) * 2 bytes = ~200 MB uncompressed
  const int n_events = 500;
  const int max_crystals = 10;
  const int ndp = 20000;

  WriteAMOREEvent(filename, n_events, max_crystals, ndp);
  ReadAMOREEvent(filename);

  return 0;
}
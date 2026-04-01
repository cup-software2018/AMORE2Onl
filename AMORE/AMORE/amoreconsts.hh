#pragma once

#include <string_view>
#include "OnlConsts/onlconsts.hh"
#include "OnlConsts/adcconsts.hh"

namespace AMORE {
    // Use constexpr for compile-time constants
    inline constexpr int kMINIMUMBCOUNT = 4096; // reading 4 [MB/try] from AMOREADC
    inline constexpr int kNCHPERADC     = 16;
    inline constexpr int kNCRYSTAL      = 18;
    inline constexpr int kNADC          = 3;

    inline constexpr int kCHUNKSIZE     = kMINIMUMBCOUNT * kKILOBYTES; 
    inline constexpr int kCHUNKNDP      = kCHUNKSIZE / 64; // # of data point in a chnunk

    inline constexpr int kRECORDLENGTH   = 30000; // record length for writing waveform
    inline constexpr int kRTRECORDLENGTH = 100000;

    // Example of using std::string_view for crystal names (C++17)
    /*
    inline constexpr std::string_view kCRYSTALNAME[] = {
        "SB28", "S35", "SS68", "SE01", "SB29", "SE02"
    };
    */
} // namespace AMORE
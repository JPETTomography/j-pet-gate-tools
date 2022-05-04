#ifndef EventAnalysis_h
#define EventAnalysis_h

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <TRandom.h>

#include "Hit.h"
#include "Common.h"

namespace event_analysis{

enum EventType {
  kUnspecified = 0,
  kTrue = 1,
  kPhantomScattered = 2,
  kDetectorScattered = 3,
  kAccidental = 4
};

enum AveragingMethod {
  kUnspecifiedAM = 0,
  kCentroidWinnerNaivelyWeighted = 1,
  kCentroidWinnerEnergyWeighted = 2,
  kCentroidWinnerEnergyWeightedFirstTime = 3,
  kEnergyWinner = 4
};

void sort_hits(std::vector<Hit> &hits, std::string key);
Hit merge_hits(const std::vector<Hit> &hits, const AveragingMethod winner);

  void RunTests();
  /// returns number of hits above noise energy threshold, number of hits above Compton energy threshold, and selected hits
  std::tuple<int, int, std::vector<Hit>> select_coincident_hits(const std::vector<Hit> &hits, double compton_energy_threshold);
  /// returns number of singles above noise energy threshold, number of singles above Compton energy threshold, and selected singles
  std::tuple<int, int, std::vector<Hit>> select_coincident_singles(const std::vector<Hit> &hits, double compton_energy_threshold);
  EventType verify_type_of_coincidence(const Hit &h1,const Hit &h2);
  void print_coincidences(const std::vector<Hit>& hits);
  void analyze_event(std::vector<Hit> &hits, bool hits_are_singles = true);

};

#endif

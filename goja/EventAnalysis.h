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
  
  // returns information if among inputed 3 hits there is exactly one prompt hit and moves the position of this hit in inputed vector to the last
  std::tuple<bool, std::vector<Hit>> get_2plus1_if_present(const std::vector<Hit>& hits, double compton_e_th_prompt);

  // Below function contains simple definition of coincidence that does not take into account Rayleigh scatterings
  EventType verify_type_of_coincidence(const Hit &h1, const Hit &h2);

  // Below function contains definition of coincidence that takes into account Rayleigh scatterings,
  // the definition is based on the snippet from the CASToR software: function ComputeKindGATEEvent from
  // https://github.com/JPETTomography/castor/blob/castor-3.1.1/src/management/gDataConversionUtilities.cc
  EventType verify_type_of_coincidence_castor(const Hit &h1, const  Hit &h2);

  void print_coincidence(const Hit& h1, const Hit& h2);
  void print_triple_coincidences(const std::vector<Hit>& hits);
  void print_coincidences(const std::vector<Hit>& hits);
  void analyze_event(std::vector<Hit> &hits, bool hits_are_singles = true, bool triple_coinc = false);

};

#endif

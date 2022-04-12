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
Hit add_hits(const std::vector<Hit> &hits, const AveragingMethod winner);

class EventAnalysis {

  std::vector<Hit> coincident_hits;
  int N;
  int N0;

public :

  EventAnalysis();

  void select_coincident_hits(const std::vector<Hit>& hits);
  void select_coincident_singles(const std::vector<Hit> &hits);
  EventType verify_type_of_coincidence(const Hit &h1,const Hit &h2) const;
  void print_coincidences();
  void analyze_event(std::vector<Hit> &hits, bool hits_are_singles = true);

};

#endif

#ifndef EventAnalysis_h
#define EventAnalysis_h

using namespace std;

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

void sort_hits(vector<Hit> &hits, string key);

class EventAnalysis {

  vector<Hit> coincident_hits;
  int N;
  int N0;

public :

  EventAnalysis();

  void select_coincident_hits(vector<Hit> &hits);
  EventType verify_type_of_coincidence(Hit &h1, Hit &h2);
  void print_coincidences();
  void analyze_event(vector<Hit> &hits);

};

#endif

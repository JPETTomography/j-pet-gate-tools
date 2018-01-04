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
#include "Constants.h"

enum EventType {
    kUnspecified = 0,
    kTrue = 1,
    kPhantomScattered = 2,
    kDetectorScattered = 3,
    kAccidental = 4
};

class EventAnalysis {

    vector<Hit> compton_hits;
    int N;
    int N0;

public :

    EventAnalysis();

    void select_compton_hits(vector<Hit> hits);
    void sort_compton_hits(string key);

    // types of coincidences: true(1), phantom-scattered(2), detector-scattered(3), acci(4)
    EventType verify_type_of_coincidence(Hit, Hit);
    void print_coincidences();
    void analyze_event(vector<Hit> hits);

};

#endif

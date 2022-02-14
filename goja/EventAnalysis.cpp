#include <map>
//#include <math.h>
#include <string>
#include <sstream>

#include "EventAnalysis.h"

//================================================================================
// MATHEMATICAL AND ADDITIONAL FUNCTIONS
//================================================================================

// struct used for sorting vectors of hits (functor)

struct EntityComp {

  string property;
  EntityComp(string property) {this->property = property;}
  bool operator()(const Hit& h1, const Hit& h2) const {
    bool val = false;
    if(property == "TIME")
      val = h1.time < h2.time;
    else if (property == "EDEP")
      val = h1.edep < h2.edep;
    else if (property == "EVENTID")
      val = h1.eventID < h2.eventID;
    return val;
  }

};

void sort_hits(vector<Hit> &hits, string key) {

  sort(hits.begin(), hits.end(), EntityComp(key));

}

Hit add_hits(const std::vector<Hit> &hits, const AveragingMethod winner = kCentroidWinnerEnergyWeightedFirstTime) {

  const unsigned int N = hits.size();
  Hit h;
  h.eventID = hits[0].eventID;
  h.volumeID = hits[0].volumeID;
  for (unsigned int i=0; i<N; i++) h.edep += hits[i].edep;
  if (winner == kCentroidWinnerNaivelyWeighted) {
    for (unsigned int i=0; i<N; i++) {
      h.time += hits[i].time;
      h.posX += hits[i].posX;
      h.posY += hits[i].posY;
      h.posZ += hits[i].posZ;
    }
    h.time /= N;
    h.posX /= N;
    h.posY /= N;
    h.posZ /= N;
  }
  else if (winner == kCentroidWinnerEnergyWeighted) {
    for (unsigned int i=0; i<N; i++) {
      h.time += hits[i].time * hits[i].edep;
      h.posX += hits[i].posX * hits[i].edep;
      h.posY += hits[i].posY * hits[i].edep;
      h.posZ += hits[i].posZ * hits[i].edep;
    }
    h.time /= h.edep;
    h.posX /= h.edep;
    h.posY /= h.edep;
    h.posZ /= h.edep;
  }
  else if (winner == kCentroidWinnerEnergyWeightedFirstTime) {
    std::vector<double> times;
    for (unsigned int i = 0; i<N; i++) times.push_back(hits[i].time);
    unsigned int min_index = std::distance(times.begin(),
                                           std::min_element(times.begin(), times.end()));
    h.time = hits[min_index].time;
    for (unsigned int i=0; i<N; i++) {
      h.posX += hits[i].posX * hits[i].edep;
      h.posY += hits[i].posY * hits[i].edep;
      h.posZ += hits[i].posZ * hits[i].edep;
    }
    h.posX /= h.edep;
    h.posY /= h.edep;
    h.posZ /= h.edep;
  }
  else if (winner == kEnergyWinner) {
    std::vector<double> energies;
    for (unsigned int i = 0; i<N; i++) energies.push_back(hits[i].edep);
    unsigned int max_index = std::distance(energies.begin(),
                                           std::max_element(energies.begin(), energies.end()));
    h.time = hits[max_index].time;
    h.posX = hits[max_index].posX;
    h.posY = hits[max_index].posY;
    h.posZ = hits[max_index].posZ;
  }
  h.sourcePosX = hits[0].sourcePosX;
  h.sourcePosY = hits[0].sourcePosY;
  h.sourcePosZ = hits[0].sourcePosZ;
  h.nPhantomCompton = hits[0].nPhantomCompton;
  h.nCrystalCompton = hits[0].nCrystalCompton;
  return h;

}


//================================================================================
// DEFINITIONS OF CLASS FUNCTIONS - BASIC ANALYSIS
//================================================================================

EventAnalysis::EventAnalysis() {

  coincident_hits.clear();
  N0 = 0;
  N = 0;

}

void EventAnalysis::select_coincident_hits(vector<Hit> &hits) {

  double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3;

  N0 = hits.size();
  for (unsigned int i=0; i<hits.size(); i++) {
    if (hits[i].edep>COMPTON_E_TH) coincident_hits.push_back(hits[i]);
  }
  N = coincident_hits.size();

}

void EventAnalysis::select_coincident_singles(const std::vector<Hit> &hits) {

  const string systemType = string(getenv("GOJA_SYSTEM_TYPE"));
  const double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3;

  map<string, vector<Hit>> singles_tmp;
  for (unsigned int i=0; i<hits.size(); i++) {
    string ID;
    if (systemType == "scanner")
      ID = std::to_string(hits[i].volumeID);
    else if (systemType == "cylindricalPET")
      ID = std::to_string(hits[i].rsectorID) + '_' + std::to_string(hits[i].layerID);
    if (singles_tmp.find(ID) == singles_tmp.end() ) {
      vector<Hit> tmp;
      tmp.push_back(hits[i]);
      singles_tmp[ID] = tmp;
    } else {
      singles_tmp[ID].push_back(hits[i]);
    }
  }

  vector<Hit> singles;
  map<string, vector<Hit>>::iterator it_tmp = singles_tmp.begin();
  while(it_tmp != singles_tmp.end()) {
    singles.push_back(add_hits(it_tmp->second, kCentroidWinnerEnergyWeightedFirstTime));
    it_tmp++;
  }

  N0 = singles.size();
  for (unsigned int i=0; i<singles.size(); i++) {
    if (singles[i].edep>COMPTON_E_TH) coincident_hits.push_back(singles[i]);
  }
  N = coincident_hits.size();

}

EventType EventAnalysis::verify_type_of_coincidence(Hit &h1, Hit &h2) {

  EventType t = kUnspecified;

  if (h1.eventID==h2.eventID) { //true, phantom-scattered and detector-scattered
    if (h1.nPhantomCompton==0 and h2.nPhantomCompton==0) {
      if (h1.nCrystalCompton==1 and h2.nCrystalCompton==1) { //true
        t = kTrue;
      }
      else { //detector-scattered
        t = kDetectorScattered;
      }
    }
    else { //phantom-scattered
      t = kPhantomScattered;
    }
  }
  else { //accidental
    t = kAccidental;
  }

  return t;

}

//================================================================================
// PRINTING
//================================================================================

void EventAnalysis::print_coincidences() {

  cout.setf(ios::fixed);

  Hit h1 = coincident_hits[0];
  Hit h2 = coincident_hits[1];

  cout.precision(2);
  cout << h1.posX/10. << "\t" << h1.posY/10. << "\t" << h1.posZ/10. << "\t";
  cout.precision(1);
  cout << h1.time << "\t";

  cout.precision(2);
  cout << h2.posX/10. << "\t" << h2.posY/10. << "\t" << h2.posZ/10. << "\t";
  cout.precision(1);
  cout << h2.time << "\t";

  cout << h1.volumeID << "\t";
  cout << h2.volumeID << "\t";

  cout.precision(2);
  cout << h1.edep << "\t";
  cout << h2.edep << "\t";

  cout << verify_type_of_coincidence(h1, h2) << "\t";

  cout.precision(2);
  cout << h1.sourcePosX/10. << "\t" << h1.sourcePosY/10. << "\t" << h1.sourcePosZ/10. << "\t";
  cout << h2.sourcePosX/10. << "\t" << h2.sourcePosY/10. << "\t" << h2.sourcePosZ/10. << endl;

}

//================================================================================
// MAIN ANALYSIS FUNCTION
//================================================================================

void EventAnalysis::analyze_event(vector<Hit> &hits, bool hits_are_singles)
{

  int MAX_N = int(atof(getenv("GOJA_MAX_N")));
  int MAX_N0 = int(atof(getenv("GOJA_MAX_N0")));

  sort_hits(hits, "TIME");

  if (hits_are_singles) {
    select_coincident_singles(hits);
  }
  else {
    select_coincident_hits(hits);
  }

  if (N==MAX_N and N0<=MAX_N0) print_coincidences();

  if (DEBUG) {
    cout.setf(ios::fixed);
    cout.precision(1);
    for (unsigned int i=0; i<hits.size(); i++) {
      double t = hits[i].time;
      if (i==0) cout << endl;
      cout.precision(1);
      cout << t << "\t";
      cout.precision(2);
      cout << hits[i].edep;
      if (i>0) {
        cout.precision(1);
        cout << "\t" << (hits[i].time-hits[i-1].time);
      }
      if (i==(hits.size()-1)) {
        cout.precision(1);
        cout << "\t" << (hits[i].time-hits[0].time) << endl;
        cout << N << "\t" << N0 << endl;
      }
      cout << endl;
    }
  }

}

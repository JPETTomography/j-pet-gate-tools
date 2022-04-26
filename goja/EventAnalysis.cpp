#include <map>
#include <string>
#include <sstream>
#include <cassert>

#include "EventAnalysis.h"
#include "Hits.h"

using namespace std;

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

Hit merge_hits(const std::vector<Hit> &hits, const AveragingMethod winner = kCentroidWinnerEnergyWeightedFirstTime) {

  const unsigned int nHits = hits.size();
  Hit h;
  h.eventID = hits[0].eventID;
  h.volumeID = hits[0].volumeID;
  for (unsigned int i=0; i<nHits; i++) h.edep += hits[i].edep;

  switch (winner) {

  case kCentroidWinnerNaivelyWeighted: {
    for (unsigned int i = 0; i < nHits; i++) {
      h.time += hits[i].time;
      h.posX += hits[i].posX;
      h.posY += hits[i].posY;
      h.posZ += hits[i].posZ;
    }
    h.time /= nHits;
    h.posX /= nHits;
    h.posY /= nHits;
    h.posZ /= nHits;
    break;
  }

  case kCentroidWinnerEnergyWeighted: {
    for (unsigned int i = 0; i < nHits; i++) {
      h.time += hits[i].time * hits[i].edep;
      h.posX += hits[i].posX * hits[i].edep;
      h.posY += hits[i].posY * hits[i].edep;
      h.posZ += hits[i].posZ * hits[i].edep;
    }
    h.time /= h.edep;
    h.posX /= h.edep;
    h.posY /= h.edep;
    h.posZ /= h.edep;
    break;
  }

  case kCentroidWinnerEnergyWeightedFirstTime: {
    std::vector<double> times;
    for (unsigned int i = 0; i < nHits; i++)
      times.push_back(hits[i].time);
    unsigned int min_index = std::distance(
        times.begin(), std::min_element(times.begin(), times.end()));
    h.time = hits[min_index].time;
    for (unsigned int i = 0; i < nHits; i++) {
      h.posX += hits[i].posX * hits[i].edep;
      h.posY += hits[i].posY * hits[i].edep;
      h.posZ += hits[i].posZ * hits[i].edep;
    }
    h.posX /= h.edep;
    h.posY /= h.edep;
    h.posZ /= h.edep;
    break;
  }

  case kEnergyWinner: {
    std::vector<double> energies;
    // for (unsigned int i = 0; i<nHits; i++) energies.push_back(hits[i].edep);
    unsigned int max_index = std::distance(
        energies.begin(), std::max_element(energies.begin(), energies.end()));
    h.time = hits[max_index].time;
    h.posX = hits[max_index].posX;
    h.posY = hits[max_index].posY;
    h.posZ = hits[max_index].posZ;
    break;
  }

  default: {
    /// should never happen
    assert(1 == 0);
    break;
  }
  }
  h.sourcePosX = hits[0].sourcePosX;
  h.sourcePosY = hits[0].sourcePosY;
  h.sourcePosZ = hits[0].sourcePosZ;
  h.nPhantomCompton = hits[0].nPhantomCompton;
  h.nCrystalCompton = hits[0].nCrystalCompton;
  return h;
}

namespace event_analysis{


std::tuple<int, int, std::vector<Hit>> select_coincident_hits(const vector<Hit> &hits, double compton_energy_threshold) 
{
  std::vector<Hit> selected_hits;

  int nb_above_noise_treshold = hits.size();
  for (unsigned int i=0; i<hits.size(); i++) {
    if (hits[i].edep> compton_energy_threshold) selected_hits.push_back(hits[i]);
  }
  int nb_above_compton_threshold = selected_hits.size();
  return std::make_tuple(nb_above_noise_treshold, nb_above_compton_threshold, selected_hits);
}

std::tuple<int, int, std::vector<Hit>>  select_coincident_singles(const std::vector<Hit> &hits, double compton_energy_threshold) {

  const string systemType = string(getenv("GOJA_SYSTEM_TYPE"));

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

  /// WK: Why it is always centroid here?
  vector<Hit> singles;
  map<string, vector<Hit>>::iterator it_tmp = singles_tmp.begin();
  while(it_tmp != singles_tmp.end()) {
    singles.push_back(merge_hits(it_tmp->second, kCentroidWinnerEnergyWeightedFirstTime));
    it_tmp++;
  }
  return select_coincident_hits(singles, compton_energy_threshold);
}

EventType verify_type_of_coincidence(const Hit &h1,const  Hit &h2) {

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

void print_coincidences(const std::vector<Hit>& hits) {

  assert(hits.size() ==2);
  cout.setf(ios::fixed);

  Hit h1 = hits[0];
  Hit h2 = hits[1];

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

void analyze_event(vector<Hit> &hits, bool hits_are_singles)
{

  double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3;
  int MAX_N = int(atof(getenv("GOJA_MAX_N")));
  int MAX_N0 = int(atof(getenv("GOJA_MAX_N0")));
  int N =0;
  int N0 =0;
  sort_hits(hits, "TIME");

  std::vector<Hit> selected_hits;
  if (hits_are_singles) {
    std::tie(N0, N, selected_hits) = select_coincident_singles(hits, COMPTON_E_TH);
  }
  else {
    std::tie(N0, N, selected_hits) = select_coincident_hits(hits, COMPTON_E_TH);
  }

  /// WK: Why N==MAX_N and N0<=MAX_N0
  if (N==MAX_N and N0<=MAX_N0) print_coincidences(selected_hits);

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
}

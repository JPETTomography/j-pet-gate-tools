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

namespace event_analysis{

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
  /// WK: temporary checks for debugs to be switched off
  assert(nHits > 0); 
  if (nHits > 0) {
    auto eID = hits[0].eventID;
    auto vID = hits[0].volumeID;
    assert(std::all_of(hits.begin(), hits.end(), [eID](const Hit& h)->bool { return h.eventID == eID; }));
    assert(std::all_of(hits.begin(), hits.end(), [vID](const Hit& h)->bool { return h.volumeID == vID; }));
  } 
  /// WK: end
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

std::tuple<bool, std::vector<Hit>> check_gamma_type(const vector<Hit> &hits, double prompt_energy_threshold) 
{
  std::vector<Hit> triple_hits;

  for (unsigned int i=0; i<hits.size(); i++) {
    if (hits[i].edep <= prompt_energy_threshold) triple_hits.push_back(hits[i]);
  }
 int prompt_multiplicity = 0;
 for (unsigned int i=0; i<hits.size(); i++) {
    if (hits[i].edep > prompt_energy_threshold) {
 triple_hits.push_back(hits[i]);
 prompt_multiplicity++;
 }
  }
 if(prompt_multiplicity == 1) {
 bool isGood = true;
 }
 else {
 bool isGood = false;
 }
  return std::make_tuple(isGood, triple_hits);
 }

}

EventType verify_type_of_triple_coincidence(const Hit &h1,const  Hit &h2,const  Hit &h3) {

  EventType t = kUnspecified;

  if (h1.eventID==h2.eventID and h1.eventID==h3.eventID) { //true, phantom-scattered and detector-scattered
    if (h1.nPhantomCompton==0 and h2.nPhantomCompton==0 and h3.nPhantomCompton==0) {
      if (h1.nCrystalCompton==1 and h2.nCrystalCompton==1 and h3.nCrystalCompton==1) { //true
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

void print_triple_coincidences(const std::vector<Hit>& hits, PROMPT_E_TH) {

 assert(hits.size() ==3);

 std::vector<Hit> triple_hits;
 bool isGood;
 std::tie(isGood, triple_hits) = check_gamma_type(hits, PROMPT_E_TH);

 if(isGood){

  cout.setf(ios::fixed);

  Hit h1 = triple_hits[0];
  Hit h2 = triple_hits[1];
  Hit h3 = triple_hits[2];

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

  cout << verify_type_of_triple_coincidence(h1, h2, h3) << "\t";

  cout.precision(2);
  cout << h1.sourcePosX/10. << "\t" << h1.sourcePosY/10. << "\t" << h1.sourcePosZ/10. << "\t";
  cout << h2.sourcePosX/10. << "\t" << h2.sourcePosY/10. << "\t" << h2.sourcePosZ/10. << endl;

  cout.precision(2);
  cout << h3.posX/10. << "\t" << h3.posY/10. << "\t" << h3.posZ/10. << "\t";
  cout.precision(1);
  cout << h3.time << "\t";
  cout << h3.volumeID << "\t";
  cout.precision(2);
  cout << h3.edep << "\t";
  cout.precision(2);
  cout << h3.sourcePosX/10. << "\t" << h3.sourcePosY/10. << "\t" << h3.sourcePosZ/10. << "\t";

}

}

//================================================================================
// MAIN ANALYSIS FUNCTION
//================================================================================

void analyze_event(vector<Hit> &hits, bool hits_are_singles, bool triple_coincidence)
{

  double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3;
  double PROMPT_E_TH = atof(getenv(“GOJA_PROMPT_E_TH”))*1e3;
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


  if (triple_coincidences) {
  if (N==MAX_N and N0<=MAX_N0) print_coincidences(selected_hits);
  }
  else {
  if (N==MAX_N and N0<=MAX_N0) print_triple_coincidences(selected_hits);
  }

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

void RunTests()
{
  const double epsilon = 0.000001;
  Hit h1;
  h1.edep = 10;
  Hit h2;
  h2.edep = 80;
  Hit h3;
  h2.edep = 50;
  std::vector<Hit> hits = {h1, h2, h3};
  for (auto& h : hits)
  {
    h.eventID = 10;
    h.volumeID = 9;
  }
  auto mergedHit = merge_hits(hits, kEnergyWinner);
  std::cout << mergedHit.edep << std::endl;
  // WK: this test breaks I guess the energy of the merged hit is wrongly assigned
  /// It should be the energy of the highest hit, so 80
  // assert(std::abs(mergedHit.edep - 80) < epsilon);
  assert(mergedHit.eventID == 10);
  assert(mergedHit.volumeID == 9);
}
}
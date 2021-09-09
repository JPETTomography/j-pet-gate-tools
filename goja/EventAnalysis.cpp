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

Hit add_hits(Hit &h1, Hit &h2, string winner="centroidWinner") {
  Hit h;
  h.eventID = h1.eventID;
  h.volumeID = h1.volumeID;
  h.trackID = h1.trackID;
  h.parentID = h1.parentID;
  h.primaryID = h1.primaryID;
  h.time = min(h1.time, h2.time);
  h.edep = h1.edep + h2.edep;
  if (winner == "centroidWinner") {
    h.posX = (h1.posX+h2.posX)/2.;
    h.posY = (h1.posY+h2.posY)/2.;
    h.posZ = (h1.posZ+h2.posZ)/2.;
  }
  else if (winner == "energyWinner") {
    if (h1.edep>=h2.edep) {
      h.posX = h1.posX;
      h.posY = h1.posY;
      h.posZ = h1.posZ;
    }
    else {
      h.posX = h2.posX;
      h.posY = h2.posY;
      h.posZ = h2.posZ;
    }
  }
  h.sourcePosX = h1.sourcePosX;
  h.sourcePosY = h1.sourcePosY;
  h.sourcePosZ = h1.sourcePosZ;
  h.nPhantomCompton = h1.nPhantomCompton;
  h.nCrystalCompton = h1.nCrystalCompton;
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

  N0 = hits.size();
  double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3;
  for (unsigned int i=0; i<hits.size(); i++) {
    if (hits[i].edep>COMPTON_E_TH) coincident_hits.push_back(hits[i]);
  }
  N = coincident_hits.size();

}

void EventAnalysis::select_coincident_singles(vector<Hit> &hits) {

  bool debug = true; // TODO: rmove

  N0 = hits.size();
  if (debug) std::cout << "N0=" << N0 << std::endl;

  if (debug) for (unsigned int i=0; i<hits.size(); i++) {
      cout << hits[i].volumeID << endl;
  }

  vector<Hit> singles;
  int volumeIDlast = hits[0].volumeID;
  vector<Hit> tmp;
  tmp.push_back(hits[0]);
  for (unsigned int i=1; i<hits.size(); i++) {
    if (hits[i].volumeID == volumeIDlast) {
      if (debug) cout << "\tvolumeIDs equal" << endl;
      tmp.push_back(hits[i]);
    }
    else {
      if (debug) cout << "\tvolumeIDs different" << endl;
      if (debug) std::cout << "tmp.size()==" << tmp.size() << std::endl;
      if (tmp.size()==1) {
        singles.push_back(tmp[0]);
      }
      else if (tmp.size()>1) {
        singles.push_back(add_hits(tmp[0], tmp[1])); // TODO: simplified: assumption that there are maximally 2 photons in single volume
      }
      if (debug) std::cout << "\tNtmp=" << tmp.size() << std::endl;
      tmp.clear();
    }
    volumeIDlast = hits[i].volumeID;
  }
  if (debug) std::cout << "Nsingles=" << singles.size() << std::endl;

  double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3;
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

void EventAnalysis::analyze_event(vector<Hit> &hits, bool singles)
{

  int MAX_N = int(atof(getenv("GOJA_MAX_N")));
  int MAX_N0 = int(atof(getenv("GOJA_MAX_N0")));

  sort_hits(hits, "TIME");

  if (singles) {
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

#define Hits_cxx
#include "Hits.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <math.h>

#include "EventAnalysis.h"
#include "Common.h"

using namespace std;

LoopResults Hits::Loop(bool singles) {
  ConfigParams params;
  params.Init(singles);

  LoopResults lr;
  if (fChain == 0)
    return lr;

  Long64_t nentries = fChain->GetEntriesFast();
  vector<Hit> hits;

  auto hitIsCompton = [](const string& procName, const string& treeName)->bool { 
    if (treeName == "Hits")
      return (procName=="Compton" or procName =="compt");
    return true;
  };

  auto hitIsProper = [](const Hit& hit, const string& treeName, int pdgEncoding, double compton_e_th_0)->bool {
      bool hit_is_proper = hit.edep > compton_e_th_0; // deposited energy is larger than the noise threshold
      if (treeName == "Hits")
          hit_is_proper = hit_is_proper and pdgEncoding==22; // gamma photon
      return hit_is_proper;
  };

  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    fChain->GetEntry(jentry);
    if (ientry < 0) break;

    Hit hit;
    hit.eventID = eventID;
    hit.volumeID = volumeID[1]+1; // numbering from 1 not from 0
    hit.rsectorID = rsectorID;
    hit.moduleID = moduleID;
    hit.submoduleID = submoduleID;
    hit.crystalID = crystalID;
    hit.layerID = layerID;
    hit.time = time*1e12; // convert to ps
    hit.edep = edep*1e3; // convert to keV
    hit.posX = posX;
    hit.posY = posY;
    hit.posZ = posZ;
    hit.sourcePosX = sourcePosX;
    hit.sourcePosY = sourcePosY;
    hit.sourcePosZ = sourcePosZ;
    hit.nPhantomCompton = nPhantomCompton;
    hit.nCrystalCompton = nCrystalCompton;
    hit.nPhantomRayleigh = nPhantomRayleigh;
    hit.nCrystalRayleigh = nCrystalRayleigh;
    std::string procName = std::string(processName);

    std::string tree_name = std::string(getenv("GOJA_TREE_NAME"));

    if(hitIsCompton(procName, tree_name)) { // the photon is scattered using Compton scattering
      lr.counter_all_compton_hits += 1;
      if (hitIsProper(hit, tree_name, PDGEncoding, params.COMPTON_E_TH_0)) {
        hits.push_back(hit);
      }
    }
  }

  if (params.EVENTS_SEPARATION_USING_TIME_WINDOW) {
    event_analysis::sort_hits(hits, "TIME");
  } else {
    if (params.EVENTS_SEPARATION_USING_IDS_OF_EVENTS) {
      event_analysis::sort_hits(hits, "EVENTID");
    } else {
      // this should never happen
      assert(1 == 0);
    }
  }

  lr.counter_compton_hits_over_the_ETH0 = hits.size();
  if (DEBUG) {
    cout.setf(ios::fixed);
    cout << "lr.counter_all_compton_hits=" << lr.counter_all_compton_hits
         << endl;
    cout << "lr.counter_compton_hits_over_the_ETH0="
         << lr.counter_compton_hits_over_the_ETH0 << endl;
  }

  Hits::FindAndDumpCoincidences(hits, params, lr);
  return lr;
}

void Hits::FindAndDumpCoincidences(const std::vector<Hit> &hits, const ConfigParams& params, LoopResults& lr)
{
    vector<Hit> event;
    int start_window_eventID = 0;
    double start_window_time = 0.;
    for (const auto &hit : hits) {
      if (DEBUG)
        cout << "hit.edep=" << hit.edep << "\thit.time=" << hit.time
             << "\tstart_window_time=" << start_window_time
             << "\thit.eventID=" << hit.eventID << endl;
      if (event.size() == 0) {
        if (hit.edep > params.COMPTON_E_TH) {
          event.push_back(
              hit); // start forming the event with the hit with edep>E_TH
          start_window_time = hit.time;
          start_window_eventID = hit.eventID;
          lr.counter_compton_hits_over_the_ETH += 1;
        }
      } else {
        bool hit_in_event = false;
        if (params.EVENTS_SEPARATION_USING_TIME_WINDOW) {
          hit_in_event = (hit.time - start_window_time) < params.TIME_WINDOW;
        } else {
          if (params.EVENTS_SEPARATION_USING_IDS_OF_EVENTS)
            hit_in_event = hit.eventID == start_window_eventID;
        }
        if (hit_in_event) { // if the current hit belongs to the current event
          event.push_back(hit); // then it is added to the current event
        } else { // if the current hit does not belong to the current event
          lr.multiplicities.push_back(event.size());
          if (event.size() >=
              2) { // if the number of hits in the current event is at least 2
            
            event_analysis::analyze_event(event,
                             params.SINGLES); // then the current event is analyzed
          }
          event.clear();                 // the current event is destroyed
          if (hit.edep > params.COMPTON_E_TH) { // start forming the event with the hit
                                         // with edep>E_TH
            event.push_back(hit);
            start_window_time = hit.time;
            start_window_eventID = hit.eventID;
            lr.counter_compton_hits_over_the_ETH += 1;
          }
        }
      }
    }
    // the last event

    if (event.size() >= 1) {
      lr.multiplicities.push_back(event.size());
      if (event.size() >=
          2) { // if the number of hits in the last event is at least 2
        event_analysis::analyze_event(event, params.SINGLES); // then the last event is analyzed
      }
      event.clear(); // the last event is destroyed
    }

    lr.real_time = hits.at(hits.size() - 1).time - hits.at(0).time;
    if (DEBUG)
      cout << "lr.counter_compton_hits_over_the_ETH="
           << lr.counter_compton_hits_over_the_ETH << endl;
}


void Hits::RunTests() 
{

  Hit hit1;
  hit1.time = 10;
  Hit hit2;
  hit2.time = 20;
  vector<Hit> hits = {hit1,hit2};
  ConfigParams params;
  params.TIME_WINDOW =5;
  params.EVENTS_SEPARATION_USING_TIME_WINDOW = 1;
  params.EVENTS_SEPARATION_USING_IDS_OF_EVENTS = 0;
  LoopResults res;
  Hits::FindAndDumpCoincidences(hits,params,res);

}

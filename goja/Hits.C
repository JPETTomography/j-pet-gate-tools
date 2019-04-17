#define Hits_cxx
#include "Hits.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "EventAnalysis.h"
#include "Common.h"

LoopResults Hits::Loop() {

	double COMPTON_E_TH_0 = atof(getenv("GOJA_COMPTON_E_TH_0"))*1e3; // [COMPTON_E_TH_0]=keV
	double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"))*1e3; // [COMPTON_E_TH]=keV
	double TIME_WINDOW = atof(getenv("GOJA_TIME_WINDOW"))*1e3; // [TIME_WINDOW]=ps
	int EVENTS_SEPARATION_USING_TIME_WINDOW, EVENTS_SEPARATION_USING_IDS_OF_EVENTS;
	int sep = int(atof(getenv("GOJA_SEP")));
	if (sep==0) EVENTS_SEPARATION_USING_TIME_WINDOW=1;
	else EVENTS_SEPARATION_USING_TIME_WINDOW=0;
	EVENTS_SEPARATION_USING_IDS_OF_EVENTS = 1-EVENTS_SEPARATION_USING_TIME_WINDOW;

	LoopResults lr = {0.,vector<int>(),0,0,0};
	if (fChain == 0)
		return lr;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	vector<Hit> hits;

	for (Long64_t jentry=0; jentry<nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		Hit hit;
		hit.eventID = eventID;
		hit.volumeID = volumeID[1]+1; // numbering from 1 not from 0
		hit.trackID = trackID;
		hit.parentID = parentID;
		hit.primaryID = primaryID;
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
		string procName = string(processName);

		bool hit_is_compton = (procName=="Compton" or procName=="compt"); // the photon is scattered using Compton scattering
		if (hit_is_compton) {
			lr.counter_all_compton_hits += 1;
			bool hit_is_proper = hit.edep>COMPTON_E_TH_0 and // deposited energy is bigger than the noise threshold
			                     nPhantomRayleigh==0 and nCrystalRayleigh==0 and // the photon is not scattered using Rayleigh scattering
			                     PDGEncoding==22; // gamma photon
			if (hit_is_proper) {
				hits.push_back(hit);
			}
		}
	}

	if (EVENTS_SEPARATION_USING_TIME_WINDOW) sort_hits(hits, "TIME");
        else if (EVENTS_SEPARATION_USING_IDS_OF_EVENTS) sort_hits(hits, "EVENTID");

	lr.counter_compton_hits_over_the_ETH0 = hits.size();
	if (DEBUG) {
		cout.setf(ios::fixed);
		cout << "lr.counter_all_compton_hits=" << lr.counter_all_compton_hits << endl;
		cout << "lr.counter_compton_hits_over_the_ETH0=" << lr.counter_compton_hits_over_the_ETH0 << endl;
	}

	vector<Hit> event;
	int start_window_eventID = 0;
	double start_window_time = 0.;
	for (unsigned int i=0; i<hits.size(); i++) {

		auto hit = hits[i];
		if(DEBUG) cout << "hit.edep=" << hit.edep << "\thit.time=" << hit.time << "\tstart_window_time=" << start_window_time << "\thit.eventID=" << hit.eventID << endl;

		if (event.size()==0) {
			if (hit.edep>COMPTON_E_TH) { // start collecting data from hit with edep>E_TH and set start of time window to its time
				event.push_back(hit); // and the new event is created using the current hit (if edep>E_TH)
				start_window_time = hit.time;
				start_window_eventID = hit.eventID;
				lr.counter_compton_hits_over_the_ETH += 1;
			}
		}
		else {
			bool hit_in_event;
			if (EVENTS_SEPARATION_USING_TIME_WINDOW) hit_in_event = (hit.time-start_window_time)<TIME_WINDOW;
			else if (EVENTS_SEPARATION_USING_IDS_OF_EVENTS) hit_in_event = hit.eventID==start_window_eventID;
			if (hit_in_event) { // if the current hit belongs to the current event
				event.push_back(hit); // then it is added to the current event
			}
			else { // if the current hit does not belong to the current event
				lr.multiplicities.push_back(event.size());
				if (event.size()>=2) { // if the number of hits in the current event is at least 2
					EventAnalysis ea;
					ea.analyze_event(event); // then the current event is analyzed
				}
				event.clear(); // the current event is destroyed
				if (hit.edep>COMPTON_E_TH) { // start collecting data from hit with edep>E_TH and set start of time window to its time
					event.push_back(hit); // and the new event is created using the current hit (if edep>E_TH)
					start_window_time = hit.time;
					start_window_eventID = hit.eventID;
					lr.counter_compton_hits_over_the_ETH += 1;
				}
			}
		}
	}

	lr.real_time = hits[hits.size()-1].time - hits[0].time;
	if(DEBUG) cout << "lr.counter_compton_hits_over_the_ETH=" << lr.counter_compton_hits_over_the_ETH << endl;
	return lr;
}

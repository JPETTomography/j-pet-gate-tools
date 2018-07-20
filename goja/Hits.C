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
#include "Constants.h"

LoopResults Hits::Loop() {

	double COMPTON_E_TH_0 = atof(getenv("GOJA_COMPTON_E_TH_0"));
	double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"));
	double TIME_WINDOW = atof(getenv("GOJA_TIME_WINDOW"))/1e9; // [TIME_WINDOW]=s
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

	int last_eventID;
	double last_time = 0.;
	double start_time = 0.;
	double start_time_window = 0.;

	vector<Hit> event;

	for (Long64_t jentry=0; jentry<nentries; jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		if (jentry==0) {
			last_eventID = eventID;
			last_time = time;
			start_time = time;
		}

		Hit hit;
		hit.eventID = eventID;
		hit.volumeID = volumeID[1]+1; // numbering from 1 not from 0
		hit.trackID = trackID;
		hit.parentID = parentID;
		hit.primaryID = primaryID;
		hit.time = time;
		hit.edep = edep;
		hit.posX = posX;
		hit.posY = posY;
		hit.posZ = posZ;
		hit.sourcePosX = sourcePosX;
		hit.sourcePosY = sourcePosY;
		hit.sourcePosZ = sourcePosZ;
		hit.nPhantomCompton = nPhantomCompton;
		hit.nCrystalCompton = nCrystalCompton;
		string procName = string(processName);

		if (procName=="Compton" or procName=="compt") {
			lr.counter_all_compton_hits += 1;
			if (hit.edep>COMPTON_E_TH_0) lr.counter_compton_hits_over_the_ETH0 += 1;
			if (hit.edep>COMPTON_E_TH) lr.counter_compton_hits_over_the_ETH += 1;
		}

		bool hit_is_proper = hit.edep>COMPTON_E_TH_0 and						// deposited energy is bigger than the noise threshold
		                     (procName=="Compton" or procName=="compt") and		// the photon is scattered using Compton scattering
		                     nPhantomRayleigh==0 and nCrystalRayleigh==0 and	// the photon is not scattered using Rayleigh scattering
		                     PDGEncoding==22;									// gamma photon

		if (hit_is_proper) {						// if the current hit is proper (see above conditions)
			if (event.size()==0) {
				if (hit.edep>COMPTON_E_TH) { // start collecting data from hit with edep>E_TH and set start of time window to its time
					event.push_back(hit);
					start_time_window = time;
				}
			}
			else {
				bool hit_in_event =
					(EVENTS_SEPARATION_USING_IDS_OF_EVENTS and last_eventID==eventID) or
					(EVENTS_SEPARATION_USING_TIME_WINDOW and fabs(start_time_window-time)<TIME_WINDOW);
				if (hit_in_event) {					// if the current hit belongs to the current event
					event.push_back(hit);			// then it is added to the current event
				}
				else {								// if the current hit does not belong to the current event
					lr.multiplicities.push_back(event.size());
					if (event.size()>=2) {			// if the number of hits in the current event is at least 2
						EventAnalysis ea;
						ea.analyze_event(event);	// then the current event is analyzed
					}
					event.clear();					// the current event is destroyed
					if (hit.edep>COMPTON_E_TH) {
						event.push_back(hit);		// and the new event is created using the current hit (if edep>E_TH)
						start_time_window = time;
					}
				}
			}
		}

		last_eventID = eventID;
		last_time = time;

	}

	lr.real_time = last_time - start_time;
	return lr;
}

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

void Hits::Loop() {

	double COMPTON_E_TH_0 = atof(getenv("GOJA_COMPTON_E_TH_0"));
	double TIME_WINDOW = atof(getenv("GOJA_TIME_WINDOW"))/1e9;
	int EVENTS_SEPARATION_USING_TIME_WINDOW, EVENTS_SEPARATION_USING_IDS_OF_EVENTS;
	int sep = int(atof(getenv("GOJA_SEP")));
	if (sep==0) EVENTS_SEPARATION_USING_TIME_WINDOW=1;
	else EVENTS_SEPARATION_USING_TIME_WINDOW=0;
	EVENTS_SEPARATION_USING_IDS_OF_EVENTS = 1-EVENTS_SEPARATION_USING_TIME_WINDOW;

	if (fChain == 0)
		return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	int last_eventID;
	double last_time;

	vector<Hit> hits;
	hits.clear();

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		if (jentry==0) {
			last_eventID = eventID;
			last_time = time;
		}

		Hit h;
		h.eventID = eventID;
		h.volumeID = volumeID[1]+1; // numbering from 1 not from 0
		h.trackID = trackID;
		h.parentID = parentID;
		h.primaryID = primaryID;
		h.time = time;
		h.processName = string(processName);
		h.edep = edep;
		h.posX = posX;
		h.posY = posY;
		h.posZ = posZ;
		h.sourcePosX = sourcePosX;
		h.sourcePosY = sourcePosY;
		h.sourcePosZ = sourcePosZ;
		h.comptVolName = comptVolName;
		h.nPhantomCompton = nPhantomCompton;
		h.nCrystalCompton = nCrystalCompton;
		h.nPhantomRayleigh = nPhantomRayleigh;
		h.nCrystalRayleigh = nCrystalRayleigh;

		bool hit_is_proper = h.edep>COMPTON_E_TH_0 and
		                     (h.processName=="Compton" or h.processName=="compt") and
		                     h.nPhantomRayleigh==0 and h.nCrystalRayleigh==0 and
		                     PDGEncoding==22;

		if ((EVENTS_SEPARATION_USING_IDS_OF_EVENTS and last_eventID==eventID) or
				(EVENTS_SEPARATION_USING_TIME_WINDOW and fabs(last_time-time)<TIME_WINDOW)) {
			if (hit_is_proper) hits.push_back(h);
		}
		else {
			if (hits.size()>=2) {
				EventAnalysis ea;
				ea.analyze_event(hits);
			}
			hits.clear();
			if (hit_is_proper) hits.push_back(h);
		}

		last_eventID = eventID;
		last_time = time;

	}
}
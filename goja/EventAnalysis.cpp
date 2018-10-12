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
		return val;
	}

};

double sign(double x) {

	if (x>0)
		return 1.;
	else if (x<0)
		return -1.;
	else
		return 0;

}

//================================================================================
// DEFINITIONS OF CLASS FUNCTIONS - BASIC ANALYSIS
//================================================================================

EventAnalysis::EventAnalysis() {

	compton_hits.clear();
	N0 = 0;
	N = 0;

}

void EventAnalysis::select_compton_hits(vector<Hit> hits) {

	N0 = hits.size();
	double COMPTON_E_TH = atof(getenv("GOJA_COMPTON_E_TH"));
	for (unsigned int i=0; i<hits.size(); i++) {
		if (hits[i].edep>COMPTON_E_TH) compton_hits.push_back(hits[i]);
	}
	N = compton_hits.size();

}

void EventAnalysis::sort_compton_hits(string key) {

	sort(compton_hits.begin(), compton_hits.end(), EntityComp(key));

}

EventType EventAnalysis::verify_type_of_coincidence(Hit h1, Hit h2) {

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

	Hit h1 = compton_hits[0];
	Hit h2 = compton_hits[1];

	cout.precision(2);
	cout << h1.posX/10. << "\t" << h1.posY/10. << "\t" << h1.posZ/10. << "\t";
	cout.precision(1);
	cout << h1.time*1e12 << "\t";

	cout.precision(2);
	cout << h2.posX/10. << "\t" << h2.posY/10. << "\t" << h2.posZ/10. << "\t";
	cout.precision(1);
	cout << h2.time*1e12 << "\t";

	cout << h1.volumeID << "\t";
	cout << h2.volumeID << "\t";

	cout.precision(2);
	cout << h1.edep*1000 << "\t";
	cout << h2.edep*1000 << "\t";

	cout << verify_type_of_coincidence(h1, h2) << "\t";

	cout.precision(2);
	cout << h1.sourcePosX/10. << "\t" << h1.sourcePosY/10. << "\t" << h1.sourcePosZ/10. << "\t";
	cout << h2.sourcePosX/10. << "\t" << h2.sourcePosY/10. << "\t" << h2.sourcePosZ/10. << endl;

}

//================================================================================
// MAIN ANALYSIS FUNCTION
//================================================================================

void EventAnalysis::analyze_event(vector< Hit > hits)
{

	int MAX_N = int(atof(getenv("GOJA_MAX_N")));
	int MAX_N0 = int(atof(getenv("GOJA_MAX_N0")));

	select_compton_hits(hits);
	sort_compton_hits("TIME");

	if (N==MAX_N and N0<=MAX_N0) print_coincidences();

}

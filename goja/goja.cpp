#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>

#include "boost/program_options.hpp"

#include "Hits.h"

#define SET_GOJA_ENV_VAR(option, variable, default_value) { \
  if (vm.count(option)) { \
    stringstream ss; \
    ss << vm[option].as<string>(); \
    setenv(variable, ss.str().c_str(), 1); \
  } \
  else { \
    setenv(variable, default_value, 1); \
  } \
  }

using namespace std;
namespace po = boost::program_options;

/**
 * Following environmental variables must be set for proper working fo the
 * GOJA analyzer (they are set indirectly using goja executable options):
 *
 * GOJA_COMPTON_E_TH    fixed energy threshold [MeV]
 * GOJA_COMPTON_E_TH_0  noise energy threshold [MeV]
 * GOJA_TIME_WINDOW     time window for a coincidence [ns]
 * GOJA_MAX_N           maximum number of events above the fixed energy
 *                      threshold in the coincidence window
 * GOJA_MAX_N0          maximum number of events above the noise energy
 *                      threshold in the coincidence window
 * GOJA_SEP             separate events using time window (0) or using
 *                      IDs of hits (1)
 * GOJA_ROOT_FILENAME   path to ROOT file
 */
int main (int argc, char* argv[]) {

  bool singles = false;

  po::options_description desc("\nGOJA (GATE Output J-PET Analyzer) help\n\nAllowed options");
  desc.add_options()

  ("help", "produce help message")

  // Forming coincidences options:
  ("eth", po::value<string>(), "fixed energy threshold [MeV] (default: 0.2 MeV)")
  ("eth0", po::value<string>(), "noise energy threshold [MeV] (default: 0.0 MeV)")
  ("tw", po::value<string>(), "time window for a coincidence [ns] (default: 3 ns)")
  ("N", po::value<string>(), "maximum number of events above the fixed energy threshold in the coincidence window (default: 2)")
  ("N0", po::value<string>(), "maximum number of events above the noise energy threshold in the coincidence window (includes N, default: 1000)")
  ("sep", po::value<string>(), "separate events using time window (arg=0) or using IDs of hits (arg=1) (default: 0)")
  ("singles", po::bool_switch(&singles), "merge hits to singles")
  ("system-type", po::value<string>(), "GATE systemType: scanner or cylindricalPET")

  // Input options:
  ("root", po::value<string>(), "file path of the single GATE *.root file,"
                                "for example --root=output.root")
  ("root-many", po::value<string>(), "file path of the base GATE *.root file "
                                     "without number and extension) and number "
                                     "of files separated using comma, "
                                     "for example, option --root-many=output/output,100 will cause "
                                     "the analysis of 100 files with names from output/output1.root "
                                     "to output/output100.root "
                                     "(base name is specified in the GATE macro using the following command: "
                                     "/gate/output/root/setFileName output/output)")

  // Output options:
  ("save-real-time-to", po::value<string>(), "save real time of the simulation to the file path "
                                             "(necessary when the input consists of many files [option --root-many] "
                                             "and some of them may be broken [for example not closed properly])")
  ("save-multiplicities-to", po::value<string>(), "save multiplicities to the file path")
  ("save-statistics-to", po::value<string>(), "save statistics to the file path (only single root file mode)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (argc == 1 or vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }

  SET_GOJA_ENV_VAR("eth", "GOJA_COMPTON_E_TH", "0.2");
  SET_GOJA_ENV_VAR("eth0", "GOJA_COMPTON_E_TH_0", "0.0");
  SET_GOJA_ENV_VAR("tw", "GOJA_TIME_WINDOW", "3");
  SET_GOJA_ENV_VAR("N", "GOJA_MAX_N", "2");
  SET_GOJA_ENV_VAR("N0", "GOJA_MAX_N0", "1000");
  SET_GOJA_ENV_VAR("sep", "GOJA_SEP", "0");
  SET_GOJA_ENV_VAR("system-type", "GOJA_SYSTEM_TYPE", "scanner");

  double real_time = 0.;
  vector<int> multiplicities;
  int counter_all_compton_hits = 0;
  int counter_compton_hits_over_the_ETH0 = 0;
  int counter_compton_hits_over_the_ETH = 0;

  if (!(vm.count("root")) and !(vm.count("root-many"))) {
    cout << "You need to add root file(s) using --root or --root-many. See help.\n";
    return 1;
  }
  else if (vm.count("root")) {
    stringstream ss;
    ss << vm["root"].as<string>();
    setenv("GOJA_ROOT_FILENAME", ss.str().c_str(), 1);
    Hits h;
    LoopResults lr = h.Loop(singles);
    real_time = lr.real_time;
    multiplicities = lr.multiplicities;
    counter_all_compton_hits = lr.counter_all_compton_hits;
    counter_compton_hits_over_the_ETH0 = lr.counter_compton_hits_over_the_ETH0;
    counter_compton_hits_over_the_ETH = lr.counter_compton_hits_over_the_ETH;
  }
  else if (vm.count("root-many")) {
    stringstream ss;
    ss << vm["root-many"].as<string>();
    string substring = "";
    vector<string> substrings;
    while (getline(ss, substring, ',')) {
      substrings.push_back(substring);
    }
    string basename = substrings[0];
    int nr_of_files = stoi(substrings[1]);

    for (int i=1; i<=nr_of_files; i++) {
      string root_filename = basename + to_string(i) + ".root";
      setenv("GOJA_ROOT_FILENAME", root_filename.c_str(), 1);
      Hits h;
      LoopResults lr = h.Loop(singles);
      real_time += lr.real_time;
      multiplicities.insert(multiplicities.end(), lr.multiplicities.begin(), lr.multiplicities.end());
    }
  }

  if (vm.count("save-real-time-to")) {
    ofstream f;
    f.open(vm["save-real-time-to"].as<string>());
    f << real_time/1e12 << endl; // convert to seconds
    f.close();
  }

  if (vm.count("save-multiplicities-to")) {
    map<unsigned int, unsigned int> multiplicities_hist;
    for (unsigned int i=0; i<multiplicities.size(); i++) {
      int m = multiplicities[i];
      if (multiplicities_hist.find(m) == multiplicities_hist.end())
        multiplicities_hist[m] = 1;
      else
        multiplicities_hist[m] += 1;
    }
    ofstream f;
    f.open(vm["save-multiplicities-to"].as<string>());
    for(auto iterator = multiplicities_hist.begin(); iterator != multiplicities_hist.end(); iterator++) {
      f << iterator->first << "\t" << iterator->second << endl;
    }
    f.close();
  }

  if (vm.count("save-statistics-to")) {
    ofstream f;
    f.open(vm["save-statistics-to"].as<string>());
    f << counter_all_compton_hits << " # all Compton hits" << endl;
    f << counter_compton_hits_over_the_ETH0 << " # compton hits with edep over the ETH0" << endl;
    f << counter_compton_hits_over_the_ETH << " # compton hits with edep over the ETH" << endl;
    f.close();
  }

  return 0;
}

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>

#include "boost/program_options.hpp"

#include "Hits.h"

using namespace std;
namespace po = boost::program_options;

int main (int argc, char* argv[]) {

	po::options_description desc("\nGOJA help\n\nAllowed options");

	desc.add_options()
    ("help", "produce help message")
	("eth", po::value<double>(), "fixed energy threshold [MeV] (default: 0.2 MeV)")
    ("eth0", po::value<double>(), "noise energy threshold [MeV] (default: 0.01 MeV)")
	("tw", po::value<double>(), "time window for a coincidence [ns] (default: 3 ns)")
	("N", po::value<double>(), "maximum number of events above the fixed energy threshold in the coincidence window (default: 2)")
	("N0", po::value<double>(), "maximum number of events above the noise energy threshold in the coincidence window (default: 3)")
	("sep", po::value<int>(), "separate events using time window (arg=0) or using IDs of hits (arg=1) (default: 0)")
    ("root", po::value<string>(), "file path of the GATE *.root file");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (argc == 1 or vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}

	if (vm.count("eth")) {
		stringstream ss;
		ss << vm["eth"].as<double>();
		setenv("GOJA_COMPTON_E_TH", ss.str().c_str(), 1);
	}
	else {
		setenv("GOJA_COMPTON_E_TH", "0.2", 1);
	}

	if (vm.count("eth0")) {
		stringstream ss;
		ss << vm["eth0"].as<double>();
		setenv("GOJA_COMPTON_E_TH_0", ss.str().c_str(), 1);
	}
	else {
		setenv("GOJA_COMPTON_E_TH_0", "0.01", 1);
	}

	if (vm.count("tw")) {
		stringstream ss;
		ss << vm["tw"].as<double>();
		setenv("GOJA_TIME_WINDOW", ss.str().c_str(), 1);
	}
	else {
		setenv("GOJA_TIME_WINDOW", "3", 1);
	}

	if (vm.count("N")) {
		stringstream ss;
		ss << vm["N"].as<double>();
		setenv("GOJA_MAX_N", ss.str().c_str(), 1);
	}
	else {
		setenv("GOJA_MAX_N", "2", 1);
	}

	if (vm.count("N0")) {
		stringstream ss;
		ss << vm["N0"].as<double>();
		setenv("GOJA_MAX_N0", ss.str().c_str(), 1);
	}
	else {
		setenv("GOJA_MAX_N0", "3", 1);
	}

	if (vm.count("sep")) {
		stringstream ss;
		ss << vm["sep"].as<int>();
		setenv("GOJA_SEP", ss.str().c_str(), 1);
	}
	else {
		setenv("GOJA_SEP", "0", 1);
	}

	if (vm.count("root")) {
		stringstream ss;
		ss << vm["root"].as<string>();
		cout << ss.str().c_str();
		setenv("GOJA_ROOT_FILENAME", ss.str().c_str(), 1);
	}
	else {
		cout << "You need to add root file using --root option. See help.\n";
		return 1;
	}

	Hits h;
	h.Loop();

	return 0;
}
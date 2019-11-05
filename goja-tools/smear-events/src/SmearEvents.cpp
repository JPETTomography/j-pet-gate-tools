/*! Executable, transforms list mode GATE data
* by adding spatial and temporal blur (see SmearEvents.h)
* author: Roman Shopa
* Roman.Shopa@ncbj.gov.pl */

#include <iostream>
#include "SmearEvents.h"
#include "ParametersReader.h"

// Readout parameters (for 20/50/100 cm - long scanner)
const Readout readoutParams = { {"PMT", { {"CRT", { 0.248695, 0.314308, 0.435673 }} } },
                                {"SI", { {"CRT", { 0.178523, 0.235196, 0.365183 }} } },
                                {"WLS", {{"FWHM", {.5}} } }
                              };

int main(int argc, char *argv[])
{
    // better to set params and constants in separate file
    ParamsReader paramsReader(smearingHitsKeys);
    paramsReader.set_defaults(defaultValues);
    paramsReader.assign_values(argc,argv);

    // info message if incorrect input
    if(!(paramsReader.is_validated())){
        std::cout << "Usage: \n" << argv[0]
                  << "\t -i input_file"
                  << " -o output_file"
                  << "\t -r <readout> (PMT, SI or WLS) \n"
                  << "\t [-l <scanner length index> "
                  << " (1/2/3 - for 20/50/100 cm, respectively)]" << std::endl;
        return 0;
    }

    // set parameters and display
    std::string inputFile(paramsReader.get_parameters()["-i"]),
                outputFile(paramsReader.get_parameters()["-o"]),
                readout(paramsReader.get_parameters()["-r"]),
                scannerLength(paramsReader.get_parameters()["-l"]);
    int JPETLength = std::stoi(std::string(scannerLength))-1;

    std::cout << "\nParameters of data transformation: \n"
              << "----------------------------------" << std::endl;

    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    std::cout << "Readout: " << readout << std::endl;
    std::cout << "Scanner length ID: "
              << std::to_string(JPETLength+1) << std::endl;

    // Data smearing section
    SmearEventsController smearController(readout,
                                          readoutParams,
                                          JPETLength);
    std::cout << "------------------\n"
              << "Progress (events):" << std::endl;
    smearController.blurData(inputFile, outputFile);

    return 0;
}

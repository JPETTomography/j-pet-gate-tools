/*! SmearEvents.h
* Tools for adding some blur to temporal and spatial (Z-coordinate) 
* parameters of events, adopted to list mode ASCII format of GATE simulations
* author: Roman Shopa
* Roman.Shopa@ncbj.gov.pl */

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>
#include <unordered_map>

#include <math.h>       /* log, pow etc */
#include <random>


#ifndef SMEAREVENTS_H
#define SMEAREVENTS_H
// C_SCIN is the effective propagation velocity
// of the signal along the strip [cm/ns]
#define C_SCIN 12.6

typedef std::unordered_map <std::string,
                            std::unordered_map <std::string, std::vector<float>>
                           > Readout;
typedef std::vector<double> Vector;
// for random generator
typedef std::normal_distribution<double> Distribution;
typedef std::default_random_engine Generator;

// Readout parameters example:
//const Readout readoutParams = { {"PMT", { {"CRT", { 0.248695, 0.314308, 0.435673 }} } },
//                          {"SI", { {"CRT", { 0.178523, 0.235196, 0.365183 }} } },
//                          {"WLS", {{"FWHM", {.5}} } }
//                        };

// A structure for static generator
struct RandomGenerator {
    Distribution distributionZ;
    Distribution distributionTime;
    Generator generator;
    // constructor
    RandomGenerator(const Vector &sigmas){
        distributionZ = Distribution(0,sigmas[0]);
        distributionTime = Distribution(0,sigmas[1]*1e3); // ns to picoseconds!
    }
};

// A class for the manipulation of input data line by line
class SmearEventsController{
public:
    // constructor and destructor
    SmearEventsController(const std::string &readout,
                          const Readout &readoutHash,
                          const int &scanner_length_ID = 1)
        : m_readout(readout),
          m_addedBlur({0.,0.,0.,0.}){
        // first estimate sigmas (standard deviations for Z & t)
        m_sigmas = m_get_sigmas(scanner_length_ID, readoutHash);
    }
    ~SmearEventsController(){}

    // main function - imports data, smears it & exports
    void blurData(const std::string &inputFileName,
                  const std::string &outputFileName){
        // stream
        std::ifstream inputFile(inputFileName);
        inputFile.clear();
        inputFile.seekg(0, std::ios::beg); // move to top
        // output file (added option to delete the content if present)
        std::ofstream outputFile;
        //std::ofstream outputFile(outputFileName);
        outputFile.open(outputFileName,
                        std::ofstream::out | std::ofstream::trunc);

        // --- The LOOP ---
        // one row from data file:
        std::string textRow;
        // set random generator
        RandomGenerator randGen(m_sigmas);
        // iterator & event counter
        int i(0);
        std::stringstream s_EventNumber;

        while(std::getline(inputFile, textRow)){
            i += 1;
            // parse into values
            Vector dataRow(m_import_text_event(textRow));
            for(int i=0; i<2; i++){
                m_addedBlur[i*2] = randGen.distributionZ(randGen.generator);
                m_addedBlur[i*2+1] = randGen.distributionTime(randGen.generator);
            }
            // blur the data
            Vector dataBlurred = m_blur_hits_and_times(dataRow,
                                                       m_addedBlur);
            std::string outRow(m_export_event_to_text(dataBlurred));
            // export to output filename
            outputFile << outRow;
            // Counter for display
            s_EventNumber.str("");
            s_EventNumber << i;
            std::cout << s_EventNumber.str() << std::flush << "\r";
        }
        outputFile.close();
        std::cout << "Events exported: " << i << std::endl;
    }

private:
    std::string m_readout;
    Vector m_sigmas,
           m_addedBlur; // 4-dimensional vector (for Z1,t1,Z2,t2)

    // Returns vector of {sigma_Z, sigma_t}
    // in (centimetres, nanoseconds) (!!!)
    Vector m_get_sigmas(const unsigned int &scanner_length_ID,
                        const Readout &readoutParams){
        // here, pointers from find() are used
        // since readoutParams is const unordered_map
        if(m_readout == "WLS") {
            auto search = readoutParams.find("WLS")->second.find("FWHM");
            double sigmaZ = search->second[0]/(2*sqrt(2*log(2)));
            return { sigmaZ, sigmaZ/C_SCIN };
        }
        else {
            auto search = readoutParams.find(m_readout)->second.find("CRT");
            double sigmaT =
                    search->second[scanner_length_ID]/(4*sqrt(log(2)));
            return { sigmaT*C_SCIN, sigmaT };
        }
    }

    // Creates event: Vector with values from text line,
    // separated by delimiter
    Vector m_import_text_event(const std::string &textRow,
                               const char* delimiter = "\t"){
        double doubleValue;
        std::string charValue;
        Vector oneRow;
        // text to stream:
        std::istringstream issRow(textRow);
        while ( getline(issRow, charValue, *delimiter) ){
           // another substream
           std::istringstream ssToken(charValue);
           ssToken >> doubleValue;
           oneRow.push_back(doubleValue);
        }
        return oneRow;
    }

    // adds random blur to Z-coordinates and times
    // according to m_sigmas
    Vector m_blur_hits_and_times(const Vector &event,
                                 const Vector &addedBlur){
        // assign as initial
        Vector blurredEvent(event);
        for(unsigned int i=0; i<2; i++){
            blurredEvent[i*4+2] = blurredEvent[i*4+2] + addedBlur[i*2];
            blurredEvent[i*4+3] = blurredEvent[i*4+3] + addedBlur[i*2+1];
        }
        return blurredEvent;
    }

    // aggregates Vector of one event into string
    std::string m_export_event_to_text(const Vector &event){
        // out buffer
        std::ostringstream outBuff;
        outBuff.precision(std::numeric_limits<double>::digits10-13); // 2 decimals
        outBuff << std::fixed;
        unsigned int size = event.size();
        for(unsigned int i=0; i<size; i++){
            outBuff << event[i];
            if(i<size-1) outBuff << "\t";
            else outBuff << "\n";
        }
        return outBuff.str();
    }
};

#endif // SMEAREVENTS_H

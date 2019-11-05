/*! ParametersReader.h
* A class for input parameters of some executable,
* transforms input string into key/value format of unordered_map
* author: Roman Shopa
* Roman.Shopa@ncbj.gov.pl */

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

#ifndef PARAMETERSREADER_H
#define PARAMETERSREADER_H

typedef std::unordered_map<std::string,
                           std::vector<std::string>> Keys;
typedef std::unordered_map<std::string, std::string> Hash;

// Define here the keys for the reader
const Keys smearingHitsKeys = {
    {"mandatory",{"-i", "-o", "-r"}},
    {"optional", {"-l"}}
};
// defaults for all keys (-l 2 denotes 50 cm - long scanner)
const Hash defaultValues = {{"-i","undefined"},
                            {"-o","undefined"},
                            {"-r","undefined"},
                            {"-l","2"}};

// A class for handling the input sequence of parameters
class ParamsReader{
public:
    // constructor and destructor
    ParamsReader(const Keys &keys)
        : m_validated(false) {
        // create a hash wih keys but no values
        for(auto &groupOfKeys: keys){
            for(auto &key: groupOfKeys.second){
                m_paramsHash.insert({key,"undefined"});
            }
        }
    }
    ~ParamsReader(){}

    // getter
    Hash get_parameters(){
        return m_paramsHash;
    }
    // whether all the parameters are set correctly
    bool is_validated(){
        return m_validated;
    }

    // parses input arguments and assigns to values where possible
    void assign_values(const int &argc,
                       char *argv[]){
        // initialize vector for arguments
        std::vector<std::string> inputParams;
        for(int i=1; i<argc; i++){
            inputParams.push_back(argv[i]);
        }
        // loop through allowed keys
        for(auto &key : get_keys(m_paramsHash)){
            // difference between pointers denotes index
            std::ptrdiff_t keyIndex =
                    std::find(inputParams.begin(),
                              inputParams.end(),
                              key) - inputParams.begin();
            // if found, assign the next value
            if(unsigned(keyIndex) < inputParams.size()){
                m_paramsHash[key] = inputParams[keyIndex+1];
            }
        }
        validate_parameters();
    }

    // sets the default values
    void set_defaults(const Hash &defaults){
        for(auto position: defaults){
            // at() is used for const Hash to avoid
            // discarding qualifiers
            m_paramsHash[position.first] = defaults.at(position.first);
        }
    }

private:
    Hash m_paramsHash;
    bool m_validated;

    // gets keys from Hash
    std::vector<std::string> get_keys(const Hash &hash){
        std::vector<std::string> output;
        for(auto &position: hash){
            output.push_back(position.first);
        }
        return output;
    }
    // gets values from Hash
    std::vector<std::string> get_values(const Hash &hash){
        std::vector<std::string> output;
        for(auto &position: hash){
            output.push_back(position.second);
        }
        return output;
    }
    // set m_validated to true if all defined
    void validate_parameters(){
        std::vector<std::string> values = get_values(m_paramsHash);
        if(std::find(values.begin(),
                     values.end(),
                     "undefined") == values.end())
        m_validated = true;
    }
};

#endif // PARAMETERSREADER_H

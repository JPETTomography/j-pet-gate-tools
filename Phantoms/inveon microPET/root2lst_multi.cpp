// Convert coincidences in a root file to Siemens Preclinical listmode data file (lst)
// This _multi version reads from multiple root files
// Geron Bindseil July 2011
// Updated Sept 2012 - Added Ability to Produce Delay Packets from Random Events
// Updated Sept 2012 - Added code to skip over any event where the block difference was 3 or less (as the scanner does)
// Updated Sept 2012 - Added code to populate the time difference bits randomly according to the empirical distribution

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iomanip>


#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TRandom2.h>

using namespace std;

void printHelpAndQuit( string const* msg )
{
cerr << *msg << std::endl;
cerr << "Usage: root2lst_multi <trues> <num_root_files> <rootfile_path> <output_file_name>" << endl;
cerr << "<trues> set to 'trues' if you want only trues, 'prompts' if you want all prompts or 'prompts_plus_delays' to include delay packets from randoms" << endl;
cerr << "<num_root_files> set the number of root files you have in the directory to process" << endl;
cerr << "<first_root_file_name> set the path to the first root file. For example: /Data/inveonPET/inveonPET_0.root" << endl;
cerr << "<output_file_name> set the output file name. For example: '/Data/inveonPET/filename.lst'" << endl;
exit( EXIT_FAILURE );
}

int determineCount( Long_t counting_bits_input )
{
switch (counting_bits_input % 8) {
case 0:
    return (0);
case 1:
    return (1);
case 2:
    return (3);
case 3:
    return (2);
case 4:
    return (6);
case 5:
    return (7);
case 6:
    return (5);
case 7:
    return (4);
}		
}

int main( int argc, char** argv )
{

if( argc < 5 )
{
string const msg = "Arguments missing!!!";
printHelpAndQuit( &msg );
}

cout << "Number of Arguments = " << argc << endl;

//Set trues only to 1 if you want (default is to take all prompts).
Int_t trues_only_flag = 0;
Int_t insert_delays_flag = 0;
if( strcmp(argv[1], "trues") == 0)
{
cout << "Only trues will be written to file."<<endl;
trues_only_flag = 1;
insert_delays_flag = 0;
}
else if (strcmp(argv[1], "prompts") == 0)
{
cout << "All prompts will be written to file." << endl;
trues_only_flag = 0;
insert_delays_flag = 0;
}
else if (strcmp(argv[1], "prompts_plus_delays") == 0)
{
cout << "All prompts will be written to file." << endl;
cout << "Randoms will be used to generate delay packets." << endl;
trues_only_flag = 0;
insert_delays_flag = 1;
}
else
{
string const msg = "Argument for 'trues' or 'prompts' or 'prompts_plus_delays' missing!!!";
printHelpAndQuit( &msg );
}

gROOT->Reset();

TStopwatch timer;
timer.Start(); 
//
//Declaration of leaves types - TTree Coincidences
//  
Int_t           comptonPhantom1;
Int_t           comptonPhantom2;
Int_t           comptonCrystal1;
Int_t           comptonCrystal2;   
Int_t           crystalID1;
Int_t           crystalID2;
Int_t           blockID1;
Int_t           blockID2;
Float_t         energy1;
Float_t         energy2;   
Int_t           eventID1;
Int_t           eventID2;
Double_t         time1;
Double_t         time2;

// Initialize Top-Level Variables
Double_t Nbr_Coinc_Prompt = 0. ;
Double_t Nbr_Coinc_Random = 0. ;
Double_t Nbr_Coinc_Scatter = 0. ;
Double_t Nbr_Coinc_NotTrues = 0.;
Double_t Nbr_Coinc_Trues = 0. ;
Double_t Nbr_Coinc_WrittenToFile = 0.;
Double_t Nbr_Coinc_BlockDiffLessThan3 = 0;

Double_t time_stamp = 0.;
Float_t time_stamp_increment = 0.0002;

Int_t ringnum1 = 0;
Int_t ringposition1 = 0;
Int_t crystalposition1 = 0;
Int_t ringnum2 = 0;
Int_t ringposition2 = 0;
Int_t crystalposition2 = 0;

Int_t ringpositiondifference = 0;

// Increments for each data packet written (includes time stamp & coincidence packets)
Long_t counting_bits = 0;

//Siemens files start from a random arbitrary time step counter value.
Long_t time_step_counter = 40266905;

// These variables store the binary values to write to the data packet.
unsigned long num = 0;
unsigned short inveonEnergy1 = 0;
unsigned short inveonEnergy2 = 0;
unsigned short PD = 1;
unsigned short TD = 0;
unsigned short count;

// Initialize the random number
TRandom2 *r2 = new TRandom2();
Double_t random_num = 0.;

// Open the output file
ofstream file(argv[4], ios::binary);

// Loop over the root files

int nfiles = atoi(argv[2]);
string basefilename = argv[3];
basefilename.erase(basefilename.length()-6);
string tempstring;
stringstream currentfilename;

//Begin Loop
for (int j=0; j<nfiles; j++) {

currentfilename.str("");
currentfilename << basefilename << j << ".root";

tempstring = currentfilename.str();

TFile *f = new TFile(tempstring.c_str(), "READ" );

if( f->IsZombie() )
{
    cerr << "Cannot open the file '" << tempstring
    << "'" << endl;
    exit( EXIT_FAILURE );
}

TTree *Coincidences = (TTree*)gDirectory->Get("Coincidences");

//   
//Set branch addresses - TTree Coincicences
//  

Coincidences->SetBranchAddress("comptonPhantom1",&comptonPhantom1);
Coincidences->SetBranchAddress("comptonPhantom2",&comptonPhantom2);
Coincidences->SetBranchAddress("crystalID1",&crystalID1);
Coincidences->SetBranchAddress("crystalID2",&crystalID2);
Coincidences->SetBranchAddress("blockID1",&blockID1);
Coincidences->SetBranchAddress("blockID2",&blockID2);
Coincidences->SetBranchAddress("eventID1",&eventID1);
Coincidences->SetBranchAddress("eventID2",&eventID2);
Coincidences->SetBranchAddress("time1",&time1);
Coincidences->SetBranchAddress("time2",&time2);

Int_t nentries = Coincidences->GetEntries();
Nbr_Coinc_Prompt += nentries;
cout<<" Starting file " << j+1 << " of " << nfiles << " with number of events: "<<  nentries <<endl;
Int_t nbytes = 0;

//
//loop on the events in the TTree Coincidences
//

for (Int_t i=0; i<nentries;i++) {
    nbytes += Coincidences->GetEntry(i);
    
    num = 0; // Reset the integer that stores the binary value to write
    
    //Evaluate Prompts, Trues and Sensitivity
    if ( eventID1 != eventID2 ) {
        Nbr_Coinc_Random++;
        
        // If we are to insert delayed coincidences from the random events do this now
        if (insert_delays_flag == 1) {
            // At this point we will insert a delay packet in front
            // of the actual random event packet with the same crystal IDs
            
            // Determine the correct count number for this packet and increment the next count
            count = determineCount(counting_bits);
            counting_bits++;
            
            // compute the ring number 0-79
            ringnum1 = (blockID1/16)*20+(crystalID1/20);
            
            // compute the ring position starting from Block 0, Crystal 0 for each ring.
            ringposition1 = ((blockID1%16)*20+(crystalID1%20) - 0 + 320)%320;
            
            // the unique crystal identifier associated with the Inveon data format.
            crystalposition1 = ringnum1*320+ringposition1; 
            
            ringnum2 = (blockID2/16)*20+(crystalID2/20); 
            ringposition2 = ((blockID2%16)*20+(crystalID2%20) - 0 + 320)%320;
            crystalposition2 = ringnum2*320+ringposition2;
            
            PD = 0; // Delay event
            random_num = r2->Rndm();
            if (random_num < 0.228) {TD = 0;}
            else if (random_num < 0.421) {TD = 1;}
            else if (random_num < 0.538) {TD = 2;}
            else if (random_num < 0.591) {TD = 3;}
            else if (random_num < 0.610) {TD = 4;}
            else if (random_num < 0.616) {TD = 5;}
            else if (random_num < 0.622) {TD = 11;}
            else if (random_num < 0.641) {TD = 12;}
            else if (random_num < 0.693) {TD = 13;}
            else if (random_num < 0.809) {TD = 14;}
            else if (random_num < 1.) {TD = 15;}
            
            // Now fill the 48-bit data packet.		
            // Start num with a clean slate
            num = 0;
            num |= 0; // 47th bit
            
            num <<= 3;
            num |= count; // 46-44
            
            num <<= 1;
            num |= 0; // 43
            
            num <<= 1;
            num |= PD; // 42
            
            num <<= 4;
            num |= TD; // 41-38
            
            num <<= 2;
            num |= inveonEnergy1; // 37-36
            
            num <<= 17;
            num |= crystalposition1; //35-19
            
            num <<= 2;
            num |= inveonEnergy2; // 18-17
            
            num <<= 17;
            num |= crystalposition2; // 16-0
            
            // Write the current packet to the file
            
            file.write((char *)(&num), 6);
            Nbr_Coinc_WrittenToFile++;
        }
    }
    
    // Reject any events that are not true unscattered coincidences
    if ( eventID1 != eventID2 || comptonPhantom1 != 0 || comptonPhantom2 != 0 ) {
        Nbr_Coinc_NotTrues++;
        if ( trues_only_flag == 1 ) {
            continue;
        }
    }
    
    // Reject any event where the azimuthal block difference is 3 or less 
    // (as the scanner actually does upon acquisition)
    if ( abs(blockID1%16 - blockID2%16) < 4 ) {
        Nbr_Coinc_BlockDiffLessThan3++;
        continue;
    }
    
    // Check if 0.2 ms have elapsed and if so write a time increment packet
    // If either time1 or time2 pass 
    // Need to continue to process all data events (so can't just i++)
    // but need to also increment the counting bits
    // Need a float that stores the most recent time value written
    
    if (time1 >= time_stamp && time2 >= time_stamp){
        //Write a timestamp packet
        
        // Determine the correct count number for this packet and increment the next count
        count = determineCount(counting_bits);
        counting_bits++;
        
        // Start num with a clean slate
        num = 0;
        num |= 0; // 47th bit
        
        num <<= 3;
        num |= count; // 46-44 (counter)
        
        num <<= 4; // 43-40
        num |= 10; // 1010 in binary
        
        num <<= 40; //39-0
        num |= time_step_counter;
        
        file.write((char *)(&num), 6);
        
        time_step_counter++;
        
        // Increment the time_stamp by time_stamp_increment
        time_stamp += time_stamp_increment;
        
        // Now continue to write the coincidence data packet after this time stamp packet.
    }
    
    // At this point, a coincidence packet is ready to be written.
    
    // Determine the correct count number for this packet and increment the next count
    count = determineCount(counting_bits);
    counting_bits++;
    
    // Now convert the Gate crystal ID to the Inveon crystal ID
    
    // compute the ring number 0-79
    ringnum1 = (blockID1/16)*20+(crystalID1/20);
    
    // compute the ring position starting from Block 0, Crystal 0 for each ring.
    ringposition1 = ((blockID1%16)*20+(crystalID1%20) - 0 + 320)%320;
    
    // the unique crystal identifier associated with the Inveon data format.
    crystalposition1 = ringnum1*320+ringposition1; 
    
    ringnum2 = (blockID2/16)*20+(crystalID2/20); 
    ringposition2 = ((blockID2%16)*20+(crystalID2%20) - 0 + 320)%320;
    crystalposition2 = ringnum2*320+ringposition2;
    
    PD = 1; // Prompt Event
    
    // Set the Time difference according to the empirical distribution:
    // 0 (0.228)  1 (0.193)  2 (0.117)  3 (0.053)  4 (0.019)  5 (0.006)  6 (0.)  7 (0.)
    // 8 (0.)  9 (0.)  10 (0.)  11 (0.006)  12 (0.019)  13 (0.052)  14 (0.116)  15 (0.191)
    random_num = r2->Rndm();
    if (random_num < 0.228) {TD = 0;}
    else if (random_num < 0.421) {TD = 1;}
    else if (random_num < 0.538) {TD = 2;}
    else if (random_num < 0.591) {TD = 3;}
    else if (random_num < 0.610) {TD = 4;}
    else if (random_num < 0.616) {TD = 5;}
    else if (random_num < 0.622) {TD = 11;}
    else if (random_num < 0.641) {TD = 12;}
    else if (random_num < 0.693) {TD = 13;}
    else if (random_num < 0.809) {TD = 14;}
    else if (random_num < 1.) {TD = 15;}
    
    // Now fill the 48-bit data packet.		
    // Start num with a clean slate
    num = 0;
    num |= 0; // 47th bit
    
    num <<= 3;
    num |= count; // 46-44
    
    num <<= 1;
    num |= 0; // 43
    
    num <<= 1;
    num |= PD; // 42
    
    num <<= 4;
    num |= TD; // 41-38
    
    num <<= 2;
    num |= inveonEnergy1; // 37-36
    
    num <<= 17;
    num |= crystalposition1; //35-19
    
    num <<= 2;
    num |= inveonEnergy2; // 18-17
    
    num <<= 17;
    num |= crystalposition2; // 16-0
    
    // Write the current packet to the file
    
    file.write((char *)(&num), 6);
    Nbr_Coinc_WrittenToFile++;
    
}
delete f;
}

file.close();
timer.Stop();		

cout << "Finished processing all files." << endl;
cout << " Time: " << timer.RealTime() << " s" << endl; 

Nbr_Coinc_Trues = Nbr_Coinc_Prompt-Nbr_Coinc_NotTrues;
Nbr_Coinc_Scatter = Nbr_Coinc_NotTrues - Nbr_Coinc_Random;

cout<<""<<endl;
cout<<""<<endl;
cout<<"#P R O M P T S     = "<<Nbr_Coinc_Prompt <<" Counts"<<endl;
cout<<"#T R U E S         = "<<Nbr_Coinc_Trues  <<" Counts"<<endl;
cout<<"#S C A T T E R S   = "<<Nbr_Coinc_Scatter  <<" Counts"<<endl;
cout<<"#R A N D O M S     = "<<Nbr_Coinc_Random <<" Counts"<<endl;
cout<<"#RAND+SCATTER      = "<<Nbr_Coinc_NotTrues <<" Counts"<<endl;
cout<<"#WRITTEN TO FILE   = "<<Nbr_Coinc_WrittenToFile<<" Counts"<<endl;
cout<<""<<endl;
cout<<""<<endl;

return 0;
}
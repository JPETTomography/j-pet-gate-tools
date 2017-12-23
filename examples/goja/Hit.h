#ifndef Hit_h
#define Hit_h

#include <string>
#include <iostream>

using namespace std;

class Hit {

public :

  int eventID;          // ID of the event
  int volumeID;         // ID of the detecting volume
  int trackID;          // ID of the track
  int parentID;         // ID of the parent particle
  int primaryID;        // ID of the primary particle
  double time;          // time from the start of the simulation [s]
  double edep;          // deposited energy [MeV]
  double posX;          // absolute x position of the hit [mm]
  double posY;          // absolute y position of the hit [mm]
  double posZ;          // absolute z position of the hit [mm]
  double localPosX;     // relative (to the center of the detecting volume) x position of the hit [mm]
  double localPosY;     // relative (to the center of the detecting volume) y position of the hit [mm]
  double localPosZ;     // relative (to the center of the detecting volume) z position of the hit [mm]
  double sourcePosX;    // absolute x position of the source of the hit [mm]
  double sourcePosY;    // absolute y position of the source of the hit [mm]
  double sourcePosZ;    // absolute z position of the source of the hit [mm]
  string processName;   // name of physical process; same processes have different
                        // names for different versions of GATE; for example Compton
                        // may be called: ‘compt’ or ‘Compton’
  int nPhantomCompton;  // number of Compton scatterings in the phantom in a track
  int nCrystalCompton;  // number of Compton scatterings in the detecting volumes in a track
  int nPhantomRayleigh; // number of Rayleigh scatterings in the phantom in a track
  int nCrystalRayleigh; // number of Rayleigh scatterings in the detecting volumes in a track

  Hit();
  void print_hit();

};

#endif

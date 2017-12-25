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
  double sourcePosX;    // absolute x position of the source of the hit [mm]
  double sourcePosY;    // absolute y position of the source of the hit [mm]
  double sourcePosZ;    // absolute z position of the source of the hit [mm]
  int nPhantomCompton;  // number of Compton scatterings in the phantom in a track
  int nCrystalCompton;  // number of Compton scatterings in the detecting volumes in a track

  Hit();
  void print_hit();

};

#endif

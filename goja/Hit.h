#ifndef Hit_h
#define Hit_h

#include <string>
#include <iostream>

struct Hit {


  int eventID = -1;          // ID of the event
  int volumeID = -1;         // volume detection in 'scanner' GATE systemType
  int rsectorID = -1;        // volume detection in 'cylindricalPET' GATE systemType
  int layerID = -1;          // volume detection in 'cylindricalPET' GATE systemType
  double time = 0.;          // time from the start of the simulation [s]
  double edep = 0.;          // deposited energy [MeV]
  double posX = 0.;          // absolute x position of the hit [mm]
  double posY = 0.;          // absolute y position of the hit [mm]
  double posZ = 0.;          // absolute z position of the hit [mm]
  double sourcePosX = 0.;    // absolute x position of the source of the hit [mm]
  double sourcePosY = 0.;    // absolute y position of the source of the hit [mm]
  double sourcePosZ = 0.;    // absolute z position of the source of the hit [mm]
  int nPhantomCompton = -1;  // number of Compton scatterings in the phantom in a track
  int nCrystalCompton = -1;  // number of Compton scatterings in the detecting volumes in a track
  int nPhantomRayleigh = -1; // number of Rayleigh scatterings in the phantom in a track
  int nCrystalRayleigh = -1; // number of Rayleigh scatterings in the detecting volumes in a track

  Hit();
  void print_hit();

};

#endif

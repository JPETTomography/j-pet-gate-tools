#ifndef Hit_h
#define Hit_h

#include <string>
#include <iostream>

using namespace std;

class Hit {

public :

  int eventID, volumeID, trackID, parentID, primaryID;
  double time, edep;
  double posX, posY, posZ;
  double localPosX, localPosY, localPosZ;
  double sourcePosX, sourcePosY, sourcePosZ;
  string processName;
  string comptVolName;
  int nPhantomCompton, nCrystalCompton;
  int nPhantomRayleigh, nCrystalRayleigh;

  Hit();
  void print_hit();

};

#endif

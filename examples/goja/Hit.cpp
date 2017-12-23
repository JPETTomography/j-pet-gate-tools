#include "Hit.h"

Hit::Hit() {
  eventID= -1;
  volumeID= -1;
  trackID = -1;
  parentID = -1;
  primaryID = -1;
  time = -1;
  edep = -1;
  processName = "";
  posX = 0.;
  posY = 0.;
  posZ = 0.;
  localPosX = 0;
  localPosY = 0;
  localPosZ = 0;
  sourcePosX = 0.;
  sourcePosY = 0.;
  sourcePosZ = 0.;
  nPhantomCompton = -1;
  nCrystalCompton = -1;
  nPhantomRayleigh = -1;
  nCrystalRayleigh = -1;
}

void Hit::print_hit() {
  cout << "eventID=" << eventID << endl;
  cout << "volumeID=" << volumeID << endl;
  cout << "trackID=" << trackID << endl;
  cout << "parentID=" << parentID << endl;
  cout << "primaryID=" << primaryID << endl;
  cout << "time=" << time << endl;
  cout << "edep=" << edep << endl;
  cout << "processName=" << processName << endl;
  cout << "hit position = [" << posX << ", " << posY << ", " << posZ << "]" << endl;
  cout << "source position = [" << sourcePosX << ", " << sourcePosY << ", " << sourcePosZ << "]" << endl;
  cout << "nPhantomCompton=" << nPhantomCompton << endl;
  cout << "nCrystalCompton=" << nCrystalCompton << endl;
  cout << "nPhantomRayleigh=" << nPhantomRayleigh << endl;
  cout << "nCrystalRayleigh=" << nCrystalRayleigh << endl;
}

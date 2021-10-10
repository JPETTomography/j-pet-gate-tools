#include "Hit.h"

Hit::Hit() {
  eventID = -1;
  volumeID = -1;
  rsectorID = -1;
  layerID = -1;
  time = -1;
  edep = -1;
  posX = 0.;
  posY = 0.;
  posZ = 0.;
  sourcePosX = 0.;
  sourcePosY = 0.;
  sourcePosZ = 0.;
  nPhantomCompton = -1;
  nCrystalCompton = -1;
}

void Hit::print_hit() {
  cout << "eventID=" << eventID << endl;
  cout << "rsectorID=" << rsectorID << endl;
  cout << "layerID=" << layerID << endl;
  cout << "volumeID=" << volumeID << endl;
  cout << "time=" << time << endl;
  cout << "edep=" << edep << endl;
  cout << "hit position = [" << posX << ", " << posY << ", " << posZ << "]" << endl;
  cout << "source position = [" << sourcePosX << ", " << sourcePosY << ", " << sourcePosZ << "]" << endl;
  cout << "nPhantomCompton=" << nPhantomCompton << endl;
  cout << "nCrystalCompton=" << nCrystalCompton << endl;
}

#include "Hit.h"
using namespace std;

Hit::Hit() {

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
  cout << "nRayleighCompton=" << nPhantomRayleigh << endl;
  cout << "nRayleighCompton=" << nCrystalRayleigh << endl;

}

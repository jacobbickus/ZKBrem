#ifndef PGA_hh
#define PGA_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"

#include "G4ParticleGun.hh"

#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

class PGA : public G4VUserPrimaryGeneratorAction
{
public:
  PGA(int beamtype = 0, int beammode = 0, double beamoffset = 0.0, int ptype = 0);
  ~PGA();
  
  G4ParticleGun *GetSource() {return TheSource;}

  void GeneratePrimaries(G4Event *anEvent);

  double SampleResonances();

private:
  G4ParticleGun *TheSource;

  TH2D *hBrems;
  //TH2D *hSample;
  TH1D *hSample;
  TH1D *hBinary;
  TH1D *hBremsWeight;

  TRandom3 Random;

  double beamo;
  int bmode;
  int beamt;

  ULong64_t genind;
};

#endif

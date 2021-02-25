// Event information class
// Based off of TrackInformation (http://geant4.slac.stanford.edu/Tips/event/1.html)
// Jayson Vavrek, MIT, 2016
// jvavrek@mit.edu

#include "eventInformation.hh"
#include "G4ios.hh"

//G4Allocator<eventInformation> anEventInformationAllocator;

eventInformation::eventInformation()
{
  weight = 0.;
  beamEnergy = 0.;
  gammaEnergy = 0.; // might be better to have an init list
  vertexX = 0.0;
  vertexY = 0.0;
  vertexZ = 0.0;
  foilTrackInit = G4ThreeVector();
  foilPathLength = 0.0;
  platedE = 0.0;

  gEnear = 0;
  gEfar = 0;
}

eventInformation::eventInformation(const G4Event* anEvent)
{
  weight = 0.;//aTrack->GetWeight();
  beamEnergy = 0.;
  gammaEnergy = 0.;
  beamAngle = 0.;
  vertexX = 0.0;
  vertexY = 0.0;
  genind = 0;
  foilTrackInit = G4ThreeVector();
  foilPathLength = 0.0;
}

eventInformation::eventInformation(const eventInformation* anEventInfo)
{
  weight = anEventInfo->GetWeight();
  beamEnergy = anEventInfo->GetBeamEnergy();
  gammaEnergy = anEventInfo->GetGammaEnergy();
  foilTrackInit = anEventInfo->GetFoilTrackInit();
  foilPathLength = anEventInfo->GetFoilPathLength();
}

eventInformation::~eventInformation(){;}

void eventInformation::SetWeight(G4double x)
{
  weight = x;
}

void eventInformation::SetBeamEnergy(G4double x)
{
	beamEnergy = x;
}

void eventInformation::SetBeamAngle(G4double x)
{
	beamAngle = x;
}

void eventInformation::SetGammaEnergy(G4double x)
{
   gammaEnergy = x;
}

void eventInformation::SetVertexX(G4double x)
{
   vertexX = x;
}

void eventInformation::SetVertexY(G4double x)
{
   vertexY = x;
}

void eventInformation::SetVertexZ(G4double x)
{
   vertexZ = x;
}

void eventInformation::SetVertexTheta(G4double x)
{
   vertexTheta = x;
}

void eventInformation::SetVertexPhi(G4double x)
{
   vertexPhi = x;
}

void eventInformation::SetFoilTrackInit(G4ThreeVector v){
  foilTrackInit = v;
}

void eventInformation::SetFoilPathLength(G4double x){
  foilPathLength = x;
}

void eventInformation::SetGeneratorIndex(ULong64_t x){
  genind = x;
}

void eventInformation::SetPlateDepartureE(G4double x){
  platedE = x;
}

void eventInformation::Print() const
{
    //G4cout 
     //<< "Original track ID " << originalTrackID 
     //<< " at " << originalPosition << G4endl;
    ;
}

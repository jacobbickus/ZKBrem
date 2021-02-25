#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include "stackingAction.hh"
#include "rootStorageManager.hh"
#include "eventInformation.hh"


stackingAction::stackingAction()
{
  beamDir = G4ThreeVector(0,0,1);
}


stackingAction::~stackingAction()
{;}


G4ClassificationOfNewTrack stackingAction::ClassifyNewTrack(const G4Track* currentTrack)
{
  // as soon as an NRF photon is generated, save it, then kill it
  const G4VProcess *cproc = currentTrack->GetCreatorProcess();
  G4VPhysicalVolume *vol = currentTrack->GetVolume();
  G4String cname, vname;

  if (cproc != NULL && vol != NULL){
    vname = vol->GetName();
    if (vname == "AlFoilLayer0" || vname == "DUFoilLayer0"){
      cname = cproc->GetProcessName();
      if (cname == "NRF"){
        G4double trackE = currentTrack->GetKineticEnergy();

        G4ThreeVector trackPosition = currentTrack->GetPosition();
        G4double trackX = trackPosition.x();
        G4double trackY = trackPosition.y();
        G4double trackZ = trackPosition.z();

        G4ThreeVector trackMomentum = currentTrack->GetMomentum();
        G4double trackThetaOut = beamDir.angle(trackMomentum);
        G4double trackPhiOut = beamDir.azimAngle(trackMomentum); // FIXME does this work?
        G4double trackpX = trackMomentum.x();
        G4double trackpY = trackMomentum.y();
        G4double trackpZ = trackMomentum.z();

        eventInformation *info = (eventInformation*) G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation();
        G4double trackThetaIn = info->GetBeamAngle();
        G4double trackEi      = info->GetBeamEnergy();
        G4double trackWeight  = info->GetWeight();

        G4bool isDU = (vname == "DUFoilLayer0");
        if (isDU) info->SetFoilTrackInit(trackPosition);

        rootStorageManager * rm = rootStorageManager::GetInstance();
        rm->FillFoilTree(trackE, trackEi, trackX, trackY, trackZ, trackThetaIn, trackThetaOut, trackPhiOut, trackpX, trackpY, trackpZ, trackWeight, isDU);
        // return fKill; // could also let live and track to the detector if we really wanted... might be worth it
      }
    }
  }

  return fUrgent;
}

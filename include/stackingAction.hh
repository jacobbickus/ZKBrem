#ifndef stackingAction_hh
#define stackingAction_hh 1

#include "G4UserStackingAction.hh"
#include "G4Track.hh"

#include "rootStorageManager.hh"
#include "runMetadata.hh"

class stackingAction : public G4UserStackingAction
{
public:
  stackingAction();
  ~stackingAction();
  
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

  void SetRootStorageManager(rootStorageManager *rSM) {theRSManager = rSM;}
  void SetMetadata(runMetadata *md) {theRunMetadata = md;}

  G4ThreeVector beamDir;

  rootStorageManager *theRSManager;
  runMetadata *theRunMetadata;
};

#endif

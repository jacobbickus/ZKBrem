#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "MPIManager.hh"

#include "runAction.hh"
#include "runMetadata.hh"
#include "rootStorageManager.hh"


runAction::runAction(G4bool arch)
  : sequentialArchitecture(arch), parallelArchitecture(!arch), nodeRank(0),
    totalEvents(0)
{;}


runAction::~runAction()
{;}


void runAction::BeginOfRunAction(const G4Run *)
{
  MPIManager *theMPIManager = MPIManager::GetInstance();
  nodeRank = theMPIManager->GetRank();

  G4double masterEvents = theMPIManager->GetMasterEvents();
  G4double slaveEvents = theMPIManager->GetSlaveEvents();

  G4cout << "\nZK ANNOUNCEMENT: # events in master = " << masterEvents
	 << " / # events in slave = "  << slaveEvents << "\n" << G4endl;

  theMPIManager->ForceBarrier("runAction::BeginOfRunAction()");
}


void runAction::EndOfRunAction(const G4Run *currentRun)
{
  MPIManager::GetInstance()->ForceBarrier("runAction::EndOfRunAction()");

  // In parallel architecture only the master should output information
  if(nodeRank == 0){

    if(parallelArchitecture){
      totalEvents = MPIManager::GetInstance()->GetInstance()->GetTotalEvents();
    }
    else
      totalEvents = currentRun->GetNumberOfEventToBeProcessed();
  }
  // rootStorageManager::GetInstance()->GetMetadata()->SetTotalEvents(totalEvents);
}

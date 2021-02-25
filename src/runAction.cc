#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#ifdef ZK_MPI_ENABLED
#include "MPIManager.hh"
#endif

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
#ifdef ZK_MPI_ENABLED
  MPIManager *theMPIManager = MPIManager::GetInstance();
  nodeRank = theMPIManager->GetRank();
  
  G4double masterEvents = theMPIManager->GetMasterEvents();
  G4double slaveEvents = theMPIManager->GetSlaveEvents();
  
  G4cout << "\nZK ANNOUNCEMENT: # events in master = " << masterEvents 
	 << " / # events in slave = "  << slaveEvents << "\n" << G4endl;
  
  theMPIManager->ForceBarrier("runAction::BeginOfRunAction()");
#endif
}


void runAction::EndOfRunAction(const G4Run *currentRun)
{
#ifdef ZK_MPI_ENABLED
  MPIManager::GetInstance()->ForceBarrier("runAction::EndOfRunAction()");
#endif

  // In parallel architecture only the master should output information
  if(nodeRank == 0){

    if(parallelArchitecture){
#ifdef ZK_MPI_ENABLED
      totalEvents = MPIManager::GetInstance()->GetInstance()->GetTotalEvents();
#endif
    }
    else
      totalEvents = currentRun->GetNumberOfEventToBeProcessed();
  }
  // rootStorageManager::GetInstance()->GetMetadata()->SetTotalEvents(totalEvents);
}

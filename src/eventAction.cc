#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"

#include <ctime>

#ifdef ZK_MPI_ENABLED
#include "MPIManager.hh"
#endif

#include "eventAction.hh"
#include "rootStorageManager.hh"

eventAction::eventAction()
  : eventInfoFreq(100000), runID(0),
    runTime(0.), prevRunTime(0.), eventsPerSec(0.), totalEventsToRun(0.), timeToFinish(0.)
{ 
  EventActionMessenger = new eventActionMessenger(this); 
}


eventAction::~eventAction()
{ delete EventActionMessenger; }


void eventAction::BeginOfEventAction(const G4Event *currentEvent)
{
  G4int event = currentEvent->GetEventID();

  if(event==0){
    G4cout << "\n\n"
	   << "\n***************************************************************************\n"
	   <<   "****  ZK STATUS: Tracking Events!  ****************************************\n***" 
	   << G4endl;

    totalEventsToRun = 
      G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  }
  else if(event % eventInfoFreq == 0){
    
    G4RunManager *runMgr = G4RunManager::GetRunManager();
    
    // Account for the fact that the process clock runs continously
    // across successive runs by subtracting the previous run time
    // from the total process time to get the current run time
    if(runMgr->GetCurrentRun()->GetRunID()!=runID){
      prevRunTime = clock()*1.0/CLOCKS_PER_SEC;
      runID++;
    }
    
    // Calculate the rate [particles tracked / s] and the estimated
    // time to completion of the present run [m,s]
    runTime = clock()*1.0/CLOCKS_PER_SEC - prevRunTime;
    eventsPerSec = event*1.0/runTime;  // [s]
    timeToFinish = (totalEventsToRun-event)/eventsPerSec; // [s]
    
    // Output the event variables in scientific notation using
    // std::stringstreams to avoid screwing up G4cout formatting
    std::stringstream eventSS;
    eventSS.precision(3);
    eventSS << std::scientific << (double)event;

    std::stringstream tEventSS;
    tEventSS.precision(3);
    tEventSS << std::scientific << totalEventsToRun;

    G4cout << "\r**  Event [" << eventSS.str() << "/" << tEventSS.str() << "]    "
	   << std::setprecision(4)
	   << "Rate [" << eventsPerSec << "]    " 
	   << std::setprecision(2)
	   << "Time2Finish [" 
	   << ((int)timeToFinish)/3600  << "h " 
	   << ((int)timeToFinish%3600)/60 << "m " 
	   << ((int)timeToFinish%3600)%60 << "s]" 
	   << std::setprecision(6) << std::flush;
  }
}


void eventAction::EndOfEventAction(const G4Event * eve)
{
  // Get the charge deposition information
  G4int cuID = G4SDManager::GetSDMpointer()->GetCollectionID("cuScore/cuChargeDep");
  G4THitsMap<G4double>* cuHC = static_cast<G4THitsMap<G4double>*>(eve->GetHCofThisEvent()->GetHC(cuID));
  G4int auID = G4SDManager::GetSDMpointer()->GetCollectionID("auScore/auChargeDep");
  G4THitsMap<G4double>* auHC = static_cast<G4THitsMap<G4double>*>(eve->GetHCofThisEvent()->GetHC(auID));
  G4int steelID = G4SDManager::GetSDMpointer()->GetCollectionID("steelScore/steelChargeDep");
  G4THitsMap<G4double>* steelHC = static_cast<G4THitsMap<G4double>*>(eve->GetHCofThisEvent()->GetHC(steelID));
  G4int junkID = G4SDManager::GetSDMpointer()->GetCollectionID("junkScore/junkChargeDep");
  G4THitsMap<G4double>* junkHC = static_cast<G4THitsMap<G4double>*>(eve->GetHCofThisEvent()->GetHC(junkID));
  
  double screwTotal = 0.0;
  for (unsigned int j=0; j<4; j++)
  {
    char buff1[512];
    sprintf(buff1,"screwScore%d/screwChargeDep%d",j,j);
    G4int sID = G4SDManager::GetSDMpointer()->GetCollectionID(buff1);
    if (sID>=0)
    {
      G4THitsMap<G4double>* sHC = static_cast<G4THitsMap<G4double>*>(eve->GetHCofThisEvent()->GetHC(sID));
      screwTotal += SumHC(sHC);
    }
  }
  
  double cuS = SumHC(cuHC);
  double auS = SumHC(auHC);
  if (junkID>=0)
    auS += SumHC(junkHC);
    
  double stS = SumHC(steelHC)+screwTotal;
  // G4cout<<SumHC(steelHC)<<" "<<screwTotal<<"\n";
  double ppS = 0;
  G4int pipeID = G4SDManager::GetSDMpointer()->GetCollectionID("pipeScore/pipeChargeDep");
  if (pipeID>=0)
  {
    G4THitsMap<G4double>* pipeHC = static_cast<G4THitsMap<G4double>*>(eve->GetHCofThisEvent()->GetHC(pipeID));
    ppS = SumHC(pipeHC);
  }
  
  //rootStorageManager::GetInstance()->FillcdepTree(cuS,auS,stS,ppS);
      
  // Reset detector pass through flags
  rootStorageManager::GetInstance()->resetHPGeFlags();
  
}

G4double eventAction::SumHC(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0;
  std::map<G4int, G4double*>::iterator it;
  for ( it = hitsMap->GetMap()->begin(); it != hitsMap->GetMap()->end(); it++) {
    sumValue += *(it->second);
  }
  return sumValue;
}

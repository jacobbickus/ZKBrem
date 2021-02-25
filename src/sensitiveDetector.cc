#include "sensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"

#include "eventInformation.hh"
#include "rootStorageManager.hh"

#include "G4Gamma.hh"
	
sensitiveDetector::sensitiveDetector(const G4String& name) 
 : G4VSensitiveDetector(name)
{}

sensitiveDetector::~sensitiveDetector() 
{}

void sensitiveDetector::Initialize(G4HCofThisEvent*)
{
  totalEne = 0.0;
  qNRF1 = false;
  qNRF2 = false;
  qNRF3 = false;

  rm = rootStorageManager::GetInstance();
	rm->resetHPGeFlags();

//  pdone = false;
//  gE = 1e10;
//  gt = 1e10;
//  gx = 1e10;
//  gy = 1e10;

}

G4bool sensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // track weight
  trackWeight = aStep->GetTrack()->GetWeight();

  // Is track the primary? --for efficiency calcs, right?
  // Entrance energy now handled in stepping action
  if (aStep->GetTrack()->GetTrackID()==1)
	{
		if (SensitiveDetectorName == "hpgeSD_1")
		{
			//if (!rm->primary[0])
      //  rm->entE[0] = aStep->GetPreStepPoint()->GetKineticEnergy();

			rm->primary[0] = true;
		}
		else if (SensitiveDetectorName == "hpgeSD_2")
		{
			//if (!rm->primary[1] && (rm->entE[1]<0))
			//	rm->entE[1] = aStep->GetPreStepPoint()->GetKineticEnergy();

			rm->primary[1] = true;
		}
		else if (SensitiveDetectorName == "hpgeSD_3")
		{
			//if (!rm->primary[2])
			//	rm->entE[2] = aStep->GetPreStepPoint()->GetKineticEnergy();

			rm->primary[2] = true;
		}
	}
 
  // NRF tagging
  const G4VProcess *proc = aStep->GetTrack()->GetCreatorProcess();
  if (proc != NULL){
    G4String procName = proc->GetProcessName();
    //G4cout << "Event: " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << G4endl;
    //G4cout << procName << G4endl;
    if (procName == "NRF")
    {
        if (SensitiveDetectorName == "hpgeSD_1") qNRF1 = true; // if any step's track comes from NRF, set the flag
        if (SensitiveDetectorName == "hpgeSD_2") qNRF2 = true; // if any step's track comes from NRF, set the flag
        if (SensitiveDetectorName == "hpgeSD_3") qNRF3 = true; // if any step's track comes from NRF, set the flag
    }
    //G4cout << "tag: " << qNRF1 << G4endl;
  }

  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep==0.) return false;
  totalEne += edep;
  
//  if (!pdone && (SensitiveDetectorName == "NearSD" || SensitiveDetectorName == "FarSD"))
//  {
//    
//    //G4cout<<SensitiveDetectorName<<"\n\n";

//    G4StepPoint * psp = aStep->GetPostStepPoint();
//    gE = psp->GetTotalEnergy();
//    gt = psp->GetMomentum().theta();
//    gx = psp->GetPosition().x();
//    gy = psp->GetPosition().y();

//    pdone = true;
//  }    

  return true;
}

void sensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  if (totalEne != 0.0)
  {
    // Account for the weights from both the initial importance sampling and the splitting.
    // There's probably a better way to do this by using a trackInformation rather than
    // an eventInformation class, but weights are multiplicative so just find the product.
    eventInformation *info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    G4double weight = info->GetWeight();
    G4double initE = info->GetBeamEnergy();
    G4double initA = info->GetBeamAngle();
    weight *= trackWeight;
    //G4cout << ">>> EOE " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " qNRF1 = " << qNRF1 << G4endl;
  	     if (SensitiveDetectorName == "hpgeSD_1") rm->FillHPGeTree1(0,totalEne, weight, qNRF1, info); // ,vertexEnergy,primary);
  	else if (SensitiveDetectorName == "hpgeSD_2") rm->FillHPGeTree1(1,totalEne, weight, qNRF2, info); //,vertexEnergy,primary);
    else if (SensitiveDetectorName == "hpgeSD_3") rm->FillHPGeTree1(2,totalEne, weight, qNRF3, info);
  	else if (SensitiveDetectorName == "lysoSD") rm->FillLYSOTree(totalEne, weight, initE, initA);
  }
  // These following fills take care of particles that pass through the HPGe
  // detectors without interaction for efficiency calculations
	else if ((rm->primary[0] || rm->indead[0]) && SensitiveDetectorName == "hpgeSD_1")
	{
    eventInformation *info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    G4double weight = info->GetWeight();
		rm->FillHPGeTree1(0, totalEne, weight, qNRF1, info); // ,vertexEnergy,primary);
	}
	else if ((rm->primary[1] || rm->indead[1]) && SensitiveDetectorName == "hpgeSD_2")
	{
    eventInformation *info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    G4double weight = info->GetWeight();
		rm->FillHPGeTree1(1, totalEne, weight, qNRF2, info); // ,vertexEnergy,primary);
	}
	else if ((rm->primary[2] || rm->indead[2]) && SensitiveDetectorName == "hpgeSD_3")
	{
    eventInformation *info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    G4double weight = info->GetWeight();
		rm->FillHPGeTree1(2, totalEne, weight, qNRF3, info); // ,vertexEnergy,primary);
	}

//  G4cout<<SensitiveDetectorName<<"\n";

//  if (gE<1e9)
//  {
//    if (SensitiveDetectorName == "NearSD") rm->FillTree(gE,gt,gx,gy,0);
//    if (SensitiveDetectorName == "FarSD") rm->FillTree(gE,gt,gx,gy,1);
//  }
}

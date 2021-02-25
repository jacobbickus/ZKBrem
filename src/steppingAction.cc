#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4SDManager.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Gamma.hh"

#include "steppingAction.hh"
#include "rootStorageManager.hh"
#include "eventInformation.hh"


steppingAction::steppingAction()
{;}


steppingAction::~steppingAction()
{;}


void steppingAction::UserSteppingAction(const G4Step *currentStep)
{

	// Check particle volume to record pass-throughs with no energy deposition for
	// certain elements
	G4StepPoint *prevPoint = currentStep->GetPreStepPoint();
	G4StepPoint *postPoint = currentStep->GetPostStepPoint();

//  if (currentStep->GetTrack()->GetTrackID()==2 && currentStep->GetTrack()->GetCreatorProcess() != NULL)
//    G4cout<<((eventInformation*)G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation())->GetGeneratorIndex()<<" "<<currentStep->GetTrack()->GetTrackID()<<" "<<currentStep->GetPreStepPoint()->GetKineticEnergy()<<" "<<currentStep->GetPostStepPoint()->GetKineticEnergy()<<" "<<currentStep->GetTrack()->GetCreatorProcess ()->GetProcessName ()<<"\n";

  //const G4VProcess *proc = postPoint->GetProcessDefinedStep();
  //if (proc != NULL){
  //  G4String pname = proc->GetProcessName();
  //  if (pname == "Rayl"){
  //    G4cout << pname
  //           << ", dE = " << (prevPoint->GetKineticEnergy() - postPoint->GetKineticEnergy())*1e6 << " eV"
  //           << ", dpx = " << (prevPoint->GetMomentum().x() - postPoint->GetMomentum().x())
  //           << ", dpy = " << (prevPoint->GetMomentum().y() - postPoint->GetMomentum().y())
  //           << ", dpz = " << (prevPoint->GetMomentum().z() - postPoint->GetMomentum().z())
  //           << G4endl;
  //  }
  //}

  // std::cout<<"\n"<<prevPoint->GetKineticEnergy()<<" "<<prevPoint->GetTotalEnergy()<<" "<<prevPoint->GetMomentum()<<"\n";

	G4VPhysicalVolume *prevPhys = prevPoint->GetPhysicalVolume();
	G4VPhysicalVolume *postPhys = postPoint->GetPhysicalVolume();
	
	if (postPhys != NULL)
	{
		G4String prevName = prevPhys->GetName();
		G4String postName = postPhys->GetName();

    rootStorageManager * rm = rootStorageManager::GetInstance();
		
		if (postName == "HPGeCrystal1")
		  rm->inHPGe[0] = true;
		if (postName == "HPGeCrystal2")
      rm->inHPGe[1] = true;
		if (postName == "HPGeCrystal3")
      rm->inHPGe[2] = true;

    eventInformation * info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());

    if (prevName=="TheHall" && postName=="FalseNear")
    {
      // G4cout<<"Filling near\n";
      //rm->FillTree(postPoint->GetTotalEnergy(),info->GetBeamEnergy(),info->GetWeight(),postPoint->GetMomentum().theta(),postPoint->GetPosition().x(),postPoint->GetPosition().y(),currentStep->GetTrack()->GetVertexPosition().x(),currentStep->GetTrack()->GetVertexPosition().y(),currentStep->GetTrack()->GetVertexPosition().z(),0);
    }
    else if (prevName=="TheHall" && postName=="FalseFar")
    {
      // G4cout<<"\n\nHERE "<<info->GetBeamEnergy()<<"\n\n";
      rootStorageManager * rm = rootStorageManager::GetInstance();
      //rm->FillTree(postPoint->GetTotalEnergy(),info->GetBeamEnergy(),info->GetWeight(),postPoint->GetMomentum().theta(),postPoint->GetPosition().x(),postPoint->GetPosition().y(),currentStep->GetTrack()->GetVertexPosition().x(),currentStep->GetTrack()->GetVertexPosition().y(),currentStep->GetTrack()->GetVertexPosition().z(),1);
      // G4cout<<"Filling far\n";
    }
    else if (prevName.contains("DUFoil") && postName=="TheHall")
    {
      
      rootStorageManager::GetInstance()->passedU1 = true;
      //G4cout<<"Plate exit "<<currentStep->GetTrack()->GetTrackID()<<"\n";
      if (currentStep->GetTrack()->GetTrackID()==1)
      {
        //G4cout<<"recording\n";
        info->SetPlateDepartureE(currentStep->GetPostStepPoint()->GetKineticEnergy());
      }
    }
    else if (postName.contains("HPGeCrystal") && !prevName.contains("HPGeCrystal"))
    {
      rootStorageManager * rm = rootStorageManager::GetInstance();
      int dind = atoi(&postName.back());
      double pE = rm->entE[dind];
      if (currentStep->GetPostStepPoint()->GetKineticEnergy()>pE)
        rm->entE[dind] = currentStep->GetPostStepPoint()->GetKineticEnergy();
    }
    else if (postName.contains("HPGe_DL") && !prevName.contains("HPGe_DL"))
    {
      rootStorageManager * rm = rootStorageManager::GetInstance();
      int dind = atoi(&postName.back());
      rm->indead[dind] = true;
      double pE = rm->entE[dind];
      if (currentStep->GetPostStepPoint()->GetKineticEnergy()>pE)
        rm->entE[dind] = currentStep->GetPostStepPoint()->GetKineticEnergy();
    }


    // energy deposition in warhead
    if (prevName.contains("TargetLayer")){
      double edep = currentStep->GetTotalEnergyDeposit();
      eventInformation *info = (eventInformation*) G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation();
      double weight = info->GetWeight();
      //G4cout << prevName << " " << edep << G4endl;
      //if (edep > 0.0) rm->FillWarheadTree(edep, weight);
    }

    // JV version of flux incident on the foil
    if (prevName == "TheHall" && postName.contains("DUFoil") && currentStep->GetTrack()->GetParticleDefinition() == G4Gamma::Definition()){
      eventInformation *info = (eventInformation*) G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation();
      //rm->FillBremIncidentTree(prevPoint->GetTotalEnergy(), info->GetWeight());
    }

    // JV computation of path length in foil
    G4ThreeVector p = postPoint->GetMomentum();
    if (prevName.contains("DUFoil") && postName == "TheHall" && p.z() < 0){
      eventInformation *info = (eventInformation*) G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation();
      // initial position in foil is only set if NRF, so check if initPos is nonzero
      G4ThreeVector posInit = info->GetFoilTrackInit();
      if (posInit.mag() > 0){
        G4ThreeVector deltaPos = postPoint->GetPosition() - posInit;
        G4double pathLength = deltaPos.mag();
        //G4cout << "pathLength [mm] = " << pathLength << " from " << info->GetFoilTrackInit() << " to " << postPoint->GetPosition() << G4endl;
        info->SetFoilPathLength(pathLength);
      }
    }

    // JV idealDetTree tallies
    if (!prevName.contains("HPGeCrystal") && postName.contains("HPGeCrystal") && p.z() < 0){
      eventInformation *info = (eventInformation*) G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation();
      rootStorageManager * rm = rootStorageManager::GetInstance();
      G4int idealIndex = atoi(&postName.back());
      rm->FillIdealDetTree(currentStep->GetTrack()->GetKineticEnergy(), info->GetWeight(), idealIndex);
    }

	} // end postPhys check

//  // Get charge accumulated (stopped) in radiator metal
//  if (currentStep->GetTrack()->GetTrackStatus() != fAlive)
//  {
//    G4cout<<"Track stopped "<<currentStep->GetTrack()->GetTrackID()<<" "<<currentStep->GetTrack()->GetDefinition()->GetParticleName()<<" "<<postPhys->GetName()<<"  "<<currentStep->GetTrack()->GetKineticEnergy()<<"\n";
//  }


   // std::cout<<"\n"<<currentStep->GetTrack()->GetTrackID()<<"\n";

   //
   // Kill conditions for speed
   //
   // Backward going in backward region
   double pz = currentStep->GetPostStepPoint()->GetMomentum().z();
   double posz = currentStep->GetPostStepPoint()->GetPosition().z();
   if (posz<-100 && pz<0)
     currentStep->GetTrack()->SetTrackStatus(fStopAndKill);

   // End point energy
   //double E = currentStep->GetPostStepPoint()->GetTotalEnergy();
   // G4string part = currentStep->GetTrack()->GetParticleDefintion()->GetParticleName();
   //if (E<2.0)
   //   currentStep->GetTrack()->SetTrackStatus(fStopAndKill);

   // Record when you cross the planes of interest
   //rootStorageManager * rm = rootStorageManager::GetInstance();

}

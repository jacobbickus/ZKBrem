// Geant4
#include "G4SDManager.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// ROOT 
#include "TChain.h"

// C++
#include <sstream>

// ZK
#include "rootStorageManager.hh"
#include "rootStorageManagerMessenger.hh"
#include "MPIManager.hh"


rootStorageManager *rootStorageManager::theRootStorageManager = 0;


rootStorageManager *rootStorageManager::GetInstance()
{ return theRootStorageManager; }


rootStorageManager::rootStorageManager(G4bool arch, int iind = 0)
  : parallelArchitecture(!arch), MPI_Rank(0), MPI_Size(1),
    ROOTFileName(""), ROOTFile(new TFile),
    tree(new TTree), hpgeTree1(new TTree), foilTree(new TTree), warheadTree(new TTree),
    bremSampleTree(new TTree), bremIncidentTree(new TTree), idealDetTree(new TTree),
    ROOTObjectsExist(false),
    gammaEnergy(0.), gammaEinit(0.), gammaEinit2(0.), wgen(0.), wgen2(0.), gammaAngle(0.), gammaX(0.), gammaY(0.), hpgeEnergy1(0.), hpgeEnergy2(0.),
    vertexEnergy1(0.), entranceEnergy1(0.),lysoEnergy(0),lysoWeight(0),lysoInitE(0),lysoInitA(0),cuCD(0),auCD(0),stCD(0),radCD(0),primary1(false), isNRF1(false),indead1(false),
    theMetadata(),find(iind), vertexX1(0), vertexY1(0), vertexZ1(0), vertexTheta1(0), vertexPhi1(0), foilPathLength1(0.), eventIndex1(0), platedE1(0),
    foilEnergy(0.), foilEnergyInit(0.), foilX(0.), foilY(0.), foilZ(0.), foilThetaIn(0.), foilThetaOut(0.), foilPhiOut(0.), foilpX(0.), foilpY(0.), foilpZ(0.), foilWeight(0.), isDU(false),
    warheadEdep(0.), warheadWeight(0.),
    bremSampleEnergy(0.), bremSampleWeight(0.), bremIncidentEnergy(0.), bremIncidentWeight(0.),
    idealDetEnergy(0.), idealDetWeight(0.), idealDetIndex(0)
{
  if(theRootStorageManager)
    G4Exception("rootStorageManager::rootStorageManager()", 
		"ZKBremRootException001", 
		FatalException, 
		"The rootStorageManager was constructed twice!");
  else
    theRootStorageManager = this;
   
  // Trackers for detector efficiency measurements 
  resetHPGeFlags();
  
  theMessenger = new rootStorageManagerMessenger(this);
}
  
rootStorageManager::~rootStorageManager()
{
  delete ROOTFile;
  delete theMessenger; 
}

void rootStorageManager::resetHPGeFlags()
{
  inHPGe.resize(10,false);
  fHPGe.resize(10,false);
  primary.resize(10,false);
  entE.resize(10,0);
  indead.resize(10,false);

  passedU1 = false;
  platedE1 = 0.0;

  for (unsigned int j=0; j<10; ++j)
  {
    indead[j] = false;
    inHPGe[j] = false;
    fHPGe[j] = false;
    primary[j] = false;
    entE[j] = -1.0;
  }
}

void rootStorageManager::CreateROOTObjects()
{
  if(ROOTObjectsExist){
    G4cout << "\nZK ANNOUNCEMENT : ROOT objects are presently are open and writeable! ZK cannot\n"
	   <<   "                  initialize a new set until the existing ROOT objects are written to\n"
	   <<   "                  to disk via the /ZK/root/write command.\n"
	   << G4endl;
    
    return;
  }
  
#ifdef ZK_MPI_ENABLED
  MPIManager *theMPIManager = MPIManager::GetInstance();
  MPI_Rank = theMPIManager->GetRank();
  MPI_Size = theMPIManager->GetSize();
#endif

  GenerateFileNames();
  
  if(ROOTFile) delete ROOTFile;
  ROOTFile = new TFile(ROOTFileName, "recreate");
  
  tree = new TTree("GammaTreeNear","TTree to hold e- bremsstrahlung gamma events at radiator");
  tree->Branch("gammaEnergy", &gammaEnergy);
  tree->Branch("gammaEinit",&gammaEinit);
  tree->Branch("weight",&wgen);
  tree->Branch("gammaAngle", &gammaAngle);
  tree->Branch("gammaX", &gammaX);
  tree->Branch("gammaY", &gammaY);
  tree->Branch("vertexX", &vx);
  tree->Branch("vertexY", &vy);
  tree->Branch("vertexZ", &vz);

  tree2 = new TTree("GammaTreeFar","TTree to hold e- bremsstrahlung gamma events at plate");
  tree2->Branch("gammaEnergy", &gammaEnergy2);
  tree2->Branch("gammaAngle", &gammaAngle2);
  tree2->Branch("gammaX", &gammaX2);
  tree2->Branch("gammaY", &gammaY2);
  tree2->Branch("vertexX", &vx);
  tree2->Branch("vertexY", &vy);
  tree2->Branch("vertexZ", &vz);
  tree2->Branch("gammaEinit",&gammaEinit2);
  tree2->Branch("weight",&wgen2);

  hpgeTree1 = new TTree("hpgeTree1","TTree to hold Edep in HPGe detector 1");
  hpgeTree1->Branch("hpgeEnergy1", &hpgeEnergy1);
  hpgeTree1->Branch("weight1", &weight1);
  hpgeTree1->Branch("vertexEnergy1", &vertexEnergy1);
  hpgeTree1->Branch("vertexX1", &vertexX1);
  hpgeTree1->Branch("vertexY1", &vertexY1);
  hpgeTree1->Branch("vertexZ1", &vertexZ1);
  hpgeTree1->Branch("ind1",&ind1);
  hpgeTree1->Branch("vertexTheta1", &vertexTheta1);
  hpgeTree1->Branch("vertexPhi1", &vertexPhi1);
  hpgeTree1->Branch("entranceEnergy1", &entranceEnergy1);
  hpgeTree1->Branch("primary1", &primary1);
  hpgeTree1->Branch("isNRF1", &isNRF1);
  hpgeTree1->Branch("foilPathLength1", &foilPathLength1);
  hpgeTree1->Branch("eventIndex1",&eventIndex1);
  hpgeTree1->Branch("plateexitE1",&platedE1);
  hpgeTree1->Branch("passedU1",&passedU1);
  hpgeTree1->Branch("inDeadLayer1",&indead1);

  lysoTree = new TTree("lysoTree","TTree to hold Edep in LYSO detector");
  lysoTree->Branch("lysoEnergy", &lysoEnergy);
  lysoTree->Branch("lysoWeight", &lysoWeight);
  lysoTree->Branch("lysoInitE", &lysoInitE);
  lysoTree->Branch("lysoInitA", &lysoInitA);
  
  cdepTree = new TTree("cdepTree","TTree to monitor charge deposition in radiator");
  cdepTree->Branch("cuCharge", &cuCD);
  cdepTree->Branch("auCharge", &auCD);
  cdepTree->Branch("steelCharge", &stCD);
  cdepTree->Branch("pipeCharge", &ppCD);
  cdepTree->Branch("radiatorCharge", &radCD);

  foilTree = new TTree("foilTree", "TTree to hold NRF photons from stackingAction");
  foilTree->Branch("foilEnergy", &foilEnergy);
  foilTree->Branch("foilEnergyInit", &foilEnergyInit);
  foilTree->Branch("foilX", &foilX);
  foilTree->Branch("foilY", &foilY);
  foilTree->Branch("foilZ", &foilZ);
  foilTree->Branch("foilThetaIn", &foilThetaIn);
  foilTree->Branch("foilThetaOut", &foilThetaOut);
  foilTree->Branch("foilPhiOut", &foilPhiOut);
  foilTree->Branch("foilpX", &foilpX);
  foilTree->Branch("foilpY", &foilpY);
  foilTree->Branch("foilpZ", &foilpZ);
  foilTree->Branch("foilWeight", &foilWeight);
  foilTree->Branch("foilIsDU", &isDU);

  warheadTree = new TTree("warheadTree", "TTree to hold energy depositions in warhead");
  warheadTree->Branch("warheadEdep", &warheadEdep);
  warheadTree->Branch("warheadWeight", &warheadWeight);

  bremSampleTree = new TTree("bremSampleTree", "TTree to track the sampled brems photons");
  bremSampleTree->Branch("bremSampleEnergy", &bremSampleEnergy);
  bremSampleTree->Branch("bremSampleWeight", &bremSampleWeight);

  bremIncidentTree = new TTree("bremIncidentTree", "TTree to track the DU-incident brems photons");
  bremIncidentTree->Branch("bremIncidentEnergy", &bremIncidentEnergy);
  bremIncidentTree->Branch("bremIncidentWeight", &bremIncidentWeight);

  idealDetTree = new TTree("idealDetTree", "TTree to tally photons just before Pb filters");
  idealDetTree->Branch("idealDetEnergy", &idealDetEnergy);
  idealDetTree->Branch("idealDetWeight", &idealDetWeight);
  idealDetTree->Branch("idealDetIndex",  &idealDetIndex);
  
  ROOTObjectsExist = true;
}


void rootStorageManager::WriteROOTObjects(G4bool EmergencyWrite)
{
  if(!ROOTObjectsExist){
    G4cout << "\nZK ANNOUNCEMENT : The ROOT objects presently do not exist and are, therefore, unwritable!\n"
           <<   "                  They should be created via the /ZK/root/init command before\n"
	   <<   "                  events are recorded and then written to disk.\n"
	   << G4endl;

    return;
  }

  if(EmergencyWrite)
    G4cout << "\nZK ANNOUNCEMENT : The ROOT objects that presently exist are being written to disk in\n"
           <<   "                   emergency fashion to avoid losing critical data before the simulation\n"
	   <<   "                   terminates. Please issue the /ZK/root/write command before exiting\n"
	   <<   "                   ZK to avoid this message.\n"
	   << G4endl;
  
  if(parallelArchitecture){
    
    tree->Write();
    tree2->Write();
    //theMetadata->Write("theMetadata");
    hpgeTree1->Write();
    lysoTree->Write();
    cdepTree->Write();
    foilTree->Write();
    warheadTree->Write();
    bremSampleTree->Write();
    bremIncidentTree->Write();
    idealDetTree->Write();
    ROOTFile->Close();

    if(MPI_Rank == 0){
      
      G4String FinalFileName = ROOTFileName.substr(0, ROOTFileName.find(".slave"));
      TFile *FinalFile = new TFile(FinalFileName.c_str(), "update");
      FinalFile->Open(FinalFileName.c_str(), "recreate");

      // merge all the data
      // using MPI seems to break the ability to write multiple TTrees using TChain::Merge
      // instead, use the ROOT utility hadd

      //TChain *gammaTreeChain = new TChain("GammaTree");

      std::vector<G4String>::iterator it;
      G4String targets = " ";
      for(it = slaveFileNames.begin(); it != slaveFileNames.end(); it++)
      {
        //gammaTreeChain->Add((*it).c_str());
        targets += (*it).c_str();
        targets += " ";
      }
      //gammaTreeChain->Merge(FinalFileName.c_str());

      G4String haddCommand = "hadd -f " + FinalFileName + targets;
      system(haddCommand.c_str());
      
      // FinalFile->Close();
      
      for(it = slaveFileNames.begin(); it != slaveFileNames.end(); it++){
        G4String RemoveSlaveFileCmd = "rm -f " + (*it);
        system(RemoveSlaveFileCmd.c_str());
      }
    }
    ROOTObjectsExist = false;
  }
  else{
    tree->Write();
    tree2->Write();
    //theMetadata->Write("theMetadata");
    hpgeTree1->Write();
    lysoTree->Write();
    cdepTree->Write();
    foilTree->Write();
    warheadTree->Write();
    bremSampleTree->Write();
    bremIncidentTree->Write();
    idealDetTree->Write();
    //ROOTFile->Close();
    
    ROOTObjectsExist = false;
  }
}  


void rootStorageManager::FillTree(G4double gE, G4double gE0, G4double w, G4double gA,G4double x,G4double y,G4double ivx, G4double ivy, G4double ivz,unsigned int id)
{
  if (!ROOTObjectsExist) return;
  // std::cout<<"\n\nCalled "<<id<<"\n\n";

  gammaEinit = gE0;
  gammaEnergy = gE;
  gammaAngle = gA;
  gammaX = x;
  gammaY = y;
  gammaEnergy2 = gE;
  gammaEinit2 = gE0;
  gammaAngle2 = gA;
  gammaX2 = x;
  gammaY2 = y;
  vx = ivx;
  vy = ivy;
  vz = ivz;
  wgen = w;
  wgen2 = w;

  if (id==0)
   tree->Fill();
  if (id==1)
  {
   tree2->Fill();
  }
}


void rootStorageManager::FillHPGeTree1(unsigned int ind, G4double gE, G4double w, G4bool qNRF, eventInformation * ei) //, G4double gvE = 0, G4bool prim = false)
{
  if (!ROOTObjectsExist) return;

  hpgeEnergy1 = gE;
  weight1 = w;
  ind1 = ind;
  entranceEnergy1 = entE[ind];
  if (ei)
  {
    vertexEnergy1 = ei->GetBeamEnergy();
    vertexX1 = ei->GetVertexX();
    vertexY1 = ei->GetVertexY();
    vertexZ1 = ei->GetVertexZ();
    vertexTheta1 = ei->GetVertexTheta();
    vertexPhi1 = ei->GetVertexPhi();
    foilPathLength1 = ei->GetFoilPathLength();
    eventIndex1 = ei->GetGeneratorIndex();
    platedE1 = ei->GetPlateDepartureE();
  }
  isNRF1 = qNRF;
  // vertexEnergy1 = gvE;
  primary1 = primary[ind];
  indead1 = indead[ind];
  hpgeTree1->Fill();
  
  fHPGe[ind] = true;
}

void rootStorageManager::FillLYSOTree(G4double gE, G4double w, G4double Ei, G4double Ai)
{
  if (!ROOTObjectsExist) return;

  lysoEnergy = gE;
  lysoWeight = w;
  lysoInitE = Ei;
  lysoInitA = Ai; 
  // vertexEnergy1 = gvE;
  // primary1 = prim;
  lysoTree->Fill();
  
  lysofilled = true;
}

void rootStorageManager::FillcdepTree(G4double cu, G4double au, G4double st, G4double pc)
{
  if (!ROOTObjectsExist) return;

  cuCD = cu;
  auCD = au;
  stCD = st;
  ppCD = pc;
  radCD = cuCD+auCD+stCD+ppCD;
  cdepTree->Fill();
}

void rootStorageManager::FillFoilTree(G4double e, G4double ei, G4double x, G4double y, G4double z, G4double thIn, G4double thOut, G4double phiOut, G4double px, G4double py, G4double pz, G4double w, G4bool qDU)
{
  if (!ROOTObjectsExist) return;

  foilEnergy = e;
  foilEnergyInit = ei;
  foilX = x;
  foilY = y;
  foilZ = z;
  foilThetaIn  = thIn;
  foilThetaOut = thOut;
  foilPhiOut = phiOut;
  foilpX = px;
  foilpY = py;
  foilpZ = pz;
  foilWeight = w;
  isDU = qDU;
  foilTree->Fill();
}

void rootStorageManager::FillWarheadTree(G4double edep, G4double w)
{
  if (!ROOTObjectsExist) return;

  warheadEdep = edep;
  warheadWeight = w;
  warheadTree->Fill();
}

void rootStorageManager::FillBremSampleTree(G4double bse, G4double bsw){
  if (!ROOTObjectsExist) return;

  bremSampleEnergy = bse;
  bremSampleWeight = bsw;
  bremSampleTree->Fill();
}

void rootStorageManager::FillBremIncidentTree(G4double bie, G4double biw){
  if (!ROOTObjectsExist) return;

  bremIncidentEnergy = bie;
  bremIncidentWeight = biw;
  bremIncidentTree->Fill();
}

void rootStorageManager::FillIdealDetTree(G4double idE, G4double idw, G4int idi){
  if (!ROOTObjectsExist) return;

  idealDetEnergy = idE;
  idealDetWeight = idw;
  idealDetIndex  = idi;
  idealDetTree->Fill();
}

void rootStorageManager::GenerateFileNames()
{
  if(ROOTObjectsExist){
    G4cout << "\nZK ANNOUNCEMENT : The ROOT objects presently exist and have already been assigned a name!\n"
	   <<   "                  To change file names, please finish event processing and write the ROOT\n"
	   <<   "                  objects to disk via the /ZK/root/write command. Then you are free\n"
	   <<   "                  to set a new file name for subsequent event storage.\n"
	   << G4endl;
    return;
  }


  // If the user has not manually set a ROOT file name then the
  // default file name uses the date/time as a unique identifier. If
  // in parallel architecture, each slave's ROOT file is given a
  // suffix of ".slaveXXX" where XXX denotes the slave's MPI rank.

  std::stringstream ss;

  if(ROOTFileName == ""){
    time_t theTime;
    struct tm *timeInfo;
    
    time(&theTime);
    timeInfo = localtime(&theTime);
    
    const int N = 100;
    char buffer[N];

    strftime(buffer,N,"ZK_%d%b%y_%H.%M",timeInfo);

    if(parallelArchitecture)
      ss << buffer << ".root.slave" << MPI_Rank;
    else
      ss << buffer << ".root";
    
    ROOTFileName = ss.str();

    if (find>0)
    {
      char b2[256];
      sprintf(b2,"altest_%d.root",find);
      G4String fname(b2);
      ROOTFileName = fname;

    }
    
  }
  else{
    if(parallelArchitecture){
      ss << ROOTFileName << ".slave" << MPI_Rank;
      ROOTFileName = ss.str();
    }
  }

  // Generate a class member vector list (on the master node only) of
  // all the slave file names. This list will be used for MPI
  // reduction before the ROOT writing process. Note that slave file
  // names are cleared at every time CreateROOTObjects() is called in
  // order to enable multiple ROOT files to be written for many runs
  // within a single ZK session
  
  if(MPI_Rank == 0){
    if(parallelArchitecture){

      slaveFileNames.clear();

      for(int rank=0; rank<MPI_Size; rank++){
	
	size_t pos = ROOTFileName.find("slave");
	if(pos != G4String::npos){
	  ss.str("");
	  ss << ROOTFileName.substr(0,pos) << "slave" << rank;
	  slaveFileNames.push_back((G4String)ss.str());
	}
	else
	  G4Exception("rootStorageManager::BeginRunCleanup()",
		      "ZKBremRootException002",
		      FatalException,
		      "ZK could not find slaves files!");
      }
    }
  }
}


void rootStorageManager::ReduceSlaveValuesToMaster()
{
#ifdef ZK_MPI_ENABLED
  
  G4cout << "\nZK ANNOUNCEMENT: Beginning the MPI reduction of data to the master!"
	 << G4endl;

  // MPIManager *theMPImanager = MPIManager::GetInstance();

  G4cout << "\nZK ANNOUNCEMENT: Finished the MPI reduction of values to the master!\n"
	 << G4endl;
#endif
}

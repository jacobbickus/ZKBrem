#ifndef rootStorageManager_hh
#define rootStorageManager_hh 1

#include "G4Event.hh"
#include "G4Run.hh"

#include "TFile.h"
#include "TTree.h"

#include <vector>

#include "rootStorageManagerMessenger.hh"
#include "eventInformation.hh"
#include "runMetadata.hh"

class rootStorageManager
{
public:
  rootStorageManager(G4bool,int);
  ~rootStorageManager();

  static rootStorageManager *GetInstance();

  void CreateROOTObjects();
  void WriteROOTObjects(G4bool emergencyWrite=false);
  void FillTree(G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, unsigned int);
  void FillHPGeTree1(unsigned int, G4double, G4double, G4bool, eventInformation * ei =  NULL); //, G4double, G4bool);
  void FillLYSOTree(G4double,G4double,G4double,G4double);
  void FillcdepTree(G4double,G4double,G4double, G4double pc = 0);
  void FillFoilTree(G4double e, G4double ei, G4double x, G4double y, G4double z, G4double thIn, G4double thOut, G4double phiOut, G4double px, G4double py, G4double pz, G4double w, G4bool qDU);
  void FillWarheadTree(G4double edep, G4double w);
  void FillBremSampleTree(G4double bse, G4double bsw);
  void FillBremIncidentTree(G4double bie, G4double biw);
  void FillIdealDetTree(G4double idE, G4double idw, G4int idi);
  void GenerateFileNames();
  void ReduceSlaveValuesToMaster();
  void resetHPGeFlags();

  void SetFileName(G4String FN) {ROOTFileName = FN;}
  G4String GetFileName() {return ROOTFileName;}

  G4bool CheckForROOTObjects() {return ROOTObjectsExist;}

  std::vector<G4bool> primary;
  G4bool primary1;
  G4bool isNRF1;
  G4double vertexEnergy1;
  unsigned int ind1;
  G4double vertexX1;
  G4double vertexY1;
  G4double vertexZ1;
  G4double vertexTheta1;
  G4double vertexPhi1;
  G4double entranceEnergy1;
  bool passedU1;
  G4double foilPathLength1;
  ULong64_t eventIndex1;
  G4double platedE1;
  std::vector<G4double> entE;
  bool indead1;
  
  std::vector<bool> inHPGe;
  std::vector<bool> fHPGe;
  std::vector<bool> indead;
  
  bool lysofilled;
  bool inlyso;

  inline runMetadata* GetMetadata() {return theMetadata;}
  inline void SetMetadata(runMetadata *md) {theMetadata = md;}
  

private:
  static rootStorageManager *theRootStorageManager;
  G4bool parallelArchitecture;
  G4int MPI_Rank, MPI_Size;

  G4String ROOTFileName;
  std::vector<G4String> slaveFileNames;
  
  rootStorageManagerMessenger *theMessenger;

  TFile * ROOTFile;
  TTree * tree;
  TTree * tree2;
  TTree * hpgeTree1;
  TTree * lysoTree;
  TTree * cdepTree;
  TTree * foilTree;
  TTree * warheadTree;
  TTree * bremSampleTree;
  TTree * bremIncidentTree;
  TTree * idealDetTree;
  G4bool ROOTObjectsExist;

  G4double wgen, wgen2, gammaEnergy, gammaEinit, gammaAngle, gammaX, gammaY, gammaEnergy2, gammaEinit2, gammaAngle2, gammaX2, gammaY2, vx, vy, vz;
  G4double lysoEnergy, lysoWeight, lysoInitE, lysoInitA;
  G4double cuCD, auCD, stCD, ppCD, radCD;
  G4double hpgeEnergy1;
  G4double hpgeEnergy2;
  G4double weight1;
  G4double weight2;
  int find;

  G4double foilEnergy, foilEnergyInit, foilX, foilY, foilZ, foilThetaIn, foilThetaOut, foilPhiOut, foilWeight, foilpX, foilpY, foilpZ;
  G4bool isDU;

  G4double warheadEdep, warheadWeight;

  G4double bremSampleEnergy, bremIncidentEnergy, bremSampleWeight, bremIncidentWeight;

  G4double idealDetEnergy, idealDetWeight;
  G4int idealDetIndex;

  runMetadata *theMetadata;
};

#endif

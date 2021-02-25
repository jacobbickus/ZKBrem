/////////////////////////////////////////////////////////////////////////////////
//
// name: ZKBrem.cc
// date: 18 Sep 14
// auth: Zach Hartwig
// mail: hartwig@psfc.mit.edu
//
// desc: 
//
/////////////////////////////////////////////////////////////////////////////////


// Geant4 classes
#include "G4RunManager.hh" 
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIcsh.hh"
#include "G4SDManager.hh"
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4UIExecutive.hh"
#include "G4ScoringManager.hh"
//#include "G4GenericBiasingPhysics.hh"

// C++ classes
#include <ctime>
#include <unistd.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <string>
using namespace std;

// ZK Classes
#include "geometryConstruction.hh"
#include "physicsList.hh"
#include "PGA.hh"
#include "steppingAction.hh"
#include "stackingAction.hh"
#include "eventAction.hh"
#include "runAction.hh"
#include "rootStorageManager.hh"
#include "runMetadata.hh"

#ifdef ZK_MPI_ENABLED
#include "MPIManager.hh"
#endif


int main(int argc, char *argv[])
{
  if(getenv("ZKBREM_TOPDIR")==NULL){
    G4cout << "\nZK ANNOUNCEMENT: The environmental variable ZKBREM_TOPDIR must be set to\n"
	   <<   "                 the top-level ZKBrem directory in the current shell session.\n"
	   <<   "                 It is recommended to source the ZKBrem setup file 'setup.sh' in\n"
	   <<   "                 your .bashrc file!\n"
	   << G4endl;
    G4Exception("ZKBrem main()", "ZKBremMainException001", FatalException, "ZK EXCEPTION: Evasive maneuvers!\n");
  }
  
  G4String ZK_TOPDIR = (G4String)std::getenv("ZKBREM_TOPDIR");
  G4String ZK_RUNTIMEDIR = ZK_TOPDIR+"/runtime";
  

  ////////////////////////////////
  // Parse command line options //
  ////////////////////////////////

  if (argc<2)
  {
    G4cout<<"ZKBrem (vis option) (macro name) (seed/label) (LYSO/LaBr3) (plate/target) (beam E) (beam dist) (particle) (false dets) (filters) (beamoffset)\n\n"
    "- ZKBrem with no options provided will print these usage instructions\n"
    "- First option (vis option) is mandatory to pass the help message, once vis\n"
    "  option is set defaults proceed as follows, must supply all preceding arguments\n"
    "- Defaults: ZKBrem (vis option) ZKBREM_TOPDIR/runtime/ZK.mac 0 0 0 0 false true 0\n"
    "- vis option: visOff (no visualization), visOn (default Geant4 visualization),\n"
    "  visQt (Qt visualization)\n"
    "- macro name: Full path to Geant4 macro, will not be used in Qt visualization\n"
    "- Seed/label: Integer used to label output file and offset the random seed\n"
    "- LYSO/LaBr3: In-beam scintillator/shielding configuration, see documentation\n"
    "  noted at instances of 'lysoconf' in src/geometryConstruction.cc\n"
    "- plate/target: Foil/target configurations, see documentation noted at instances\n"
    "  of 'plateconf' in src/geometryConstruction.cc, uses bitwise test for the\n"
    "  Sept. 2017 configurations\n"
    "- beam E: Beam energy parameter, 0 = Uniform energy between 1.9 and 2.3 MeV\n"
    "  (range can be changed in src/PGA.cc), 1 = Sampling of U-238/Al-27 resonances\n"
    "  2 = Sept. 2017 endpoint energy fixed (2.521 MeV)\n"
    "- beam dist: Beam illumination mode, 0 = collimated bremsstrahlung beam,\n"
    "  1 = illumination of HPGe detectors originating from the foils (in their Sept.\n"
    "  2017 configuration)\n"
    "- particle: Beam particle type, 0 = gamma, 1 = geantino, 2 = electron\n"
    "- false dets: boolean for inclusion of false air detector planes immediately\n"
    "  after the radiator and before the foil used for recording the bremsstrahlung\n"
    "  production, otherwise produces masssive output files\n"
    "- filters: boolean for inclusion of lead filters between HPGe detectors and the\n"
    "  foil(s), useful for simulations of intrinsic detector efficiencies\n"
    "- beamoffset: Parameter passed to the PGA for utility, currently unused\n\n";

    // G4cout<<"\n\nZKBrem (vis option) (macro name) (seed/label) (LYSO/LaBr3) (plate/target) (beam E) (beam dist) (false dets) (filters)\n\n";

    

    return 1;
  }

  // Ability to set visualization: 'on' or 'off'
  G4bool visualization = false;
  G4bool visQt = false;
  if(argc>1){
    std::string arg1 = argv[1];
    if(arg1=="visOn")
      visualization = true;
    if(arg1=="visQt"){
      visualization = true;
      visQt = true;
    }
  }
  
  // The user may specify an alternate macro file for MPI mode;
  // otherwise, the default will be used.
  G4String MPIMacroName = ZK_RUNTIMEDIR+"/ZK.mac";
  if(argc>2)
    MPIMacroName.assign(argv[2]);

  int manualp = 0;
  if(argc>3)
    manualp = atoi(argv[3]);
    
  int plateconf = 0;
  if (argc>5)
    plateconf = atoi(argv[5]);
    
  if (plateconf<0)
  {
    G4cout<<"\n\nInvalid plate configuration option flag! ABORT!\n\n";
    exit(-1);
  }
    
  int lysoconf = 0;
  if (argc>4)
    lysoconf = atoi(argv[4]);
    
  if (lysoconf>5 || lysoconf<0)
  {
    G4cout<<"\n\nInvalid LYSO configuration option flag! ABORT!\n\n";
    exit(-2);
  }

  int beamtype = 0;
  if (argc>6)
    beamtype = atoi(argv[6]);

  int beammode = 0;
  if (argc>7)
    beammode = atoi(argv[7]);

  int particle = 0;
  if (argc>8)
    particle = atoi(argv[8]);

  bool incff = false;
  if (argc>9)
    incff = (bool)atoi(argv[9]);

  bool incfilt = true;
  if (argc>10)
    incfilt = (bool)atoi(argv[10]);

  double beamoffset = 0.0;
  if (argc>11)
    beamoffset = atof(argv[11]);

  // Select whether ZK architecture is sequential or parallel based on
  // the binary name (sequential=="ZK"), parallel=="ZK_MPI")
  
G4bool sequentialBuild = true;
  std::string arg0 = argv[0];
  if(arg0=="ZKBrem_MPI")
    sequentialBuild = false;
  
  
  //////////////////////////////////////////////
  // Initialize mandatory/user Geant4 classes //
  //////////////////////////////////////////////

  // Create the theRunManager to handle program flow
  G4RunManager *theRunManager = new G4RunManager;

    // Initialize MPI before all the user classes
#ifdef ZK_MPI_ENABLED
  const G4int argcMPI = 2;
  char *argvMPI[argcMPI];
  argvMPI[0] = argv[0]; // binary name
  argvMPI[1] = (char *)"/tmp/ZKSlave"; // slave file base name

  MPIManager *theMPIManager = new MPIManager(argcMPI,argvMPI);
#endif
  
  // Assign the mandatory user-derived classes to the run manager and
  // initialize it before creation of the user actions so that it can
  // be subsequently accessed from the constructors of the user action
  // classes. Note that PGA needs sequential/parallel in order to
  // determine if the GPS macros should be built (only in seq. arch.)
  theRunManager->SetUserInitialization(new geometryConstruction(lysoconf,plateconf,incff,incfilt));
  physicsList *thePhysics = new physicsList();
  //G4GenericBiasingPhysics *theBiasingPhysics = new G4GenericBiasingPhysics();
  // theBiasingPhysics->NonPhysicsBias("gamma");
  //thePhysics->RegisterPhysics(theBiasingPhysics);
  theRunManager->SetUserInitialization(thePhysics);

  theRunManager->SetUserAction(new PGA(beamtype,beammode,beamoffset,particle));
  theRunManager->Initialize();
  
  // Create the user action classess and assign to the run manager
  steppingAction *SteppingAction = new steppingAction();
  theRunManager->SetUserAction(SteppingAction);

  stackingAction *StackingAction = new stackingAction();
  theRunManager->SetUserAction(StackingAction);

  eventAction *EventAction = new eventAction();
  theRunManager->SetUserAction(EventAction);

  runAction *RunAction = new runAction(sequentialBuild);
  theRunManager->SetUserAction(RunAction);


  ////////////////////////////////////
  // Initialize the scoring manager //
  ////////////////////////////////////

  G4ScoringManager* theScoringManager = G4ScoringManager::GetScoringManager();
  
  
  /////////////////////////////////////////////////////////////////
  // Initialize the user-interface manager with defaults/aliases //
  /////////////////////////////////////////////////////////////////

  // Get the pointer to the U(ser) I(nterface) manager.
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  {
    G4String alias = "/control/alias";
    G4String exe = "/control/execute {run}/";

    // Set the alias to contain absolute path to runtime files
    UImanager->ApplyCommand(alias + " run " + ZK_RUNTIMEDIR);
    UImanager->ApplyCommand(alias + " mac " + exe + "ZK.mac");
    UImanager->ApplyCommand(alias + " ebeam " + exe + "electronBeamVertex.gps");
    UImanager->ApplyCommand(alias + " score " + exe + "score.mac");

    // Set useful aliases
    UImanager->ApplyCommand("/control/alias init ZK/root/init");
    UImanager->ApplyCommand("/control/alias write ZK/root/write");

    // Set some basic Geant4-specific verbosity defaults
    UImanager->ApplyCommand("/run/verbose 2");
    UImanager->ApplyCommand("/event/verbose 0");
    UImanager->ApplyCommand("/hits/verbose 0");
    UImanager->ApplyCommand("/tracking/verbose 0");
    UImanager->ApplyCommand("/control/verbose 2");
    UImanager->ApplyCommand("/gps/verbose 0");

    // Set the default particle source
    UImanager->ApplyCommand("{ebeam}");
    UImanager->ApplyCommand("/ZK/root/init");
  }  


  /////////////////////////////////
  // Initialize the ROOT manager //
  /////////////////////////////////

  rootStorageManager *theRSManager = new rootStorageManager(sequentialBuild,manualp);
  StackingAction->SetRootStorageManager(theRSManager);

  // and the metadata instance
  //runMetadata *theMetadata = new runMetadata();
  //theRSManager->SetMetadata(theMetadata);

  
  ///////////////////////////////////////
  // Current ZK instance is sequential //
  ///////////////////////////////////////

  if(sequentialBuild){
    
    // ZKBrem presently always randomizes its RNG seed; changing
    // conditional to false will start all simulations with same seed
    if(true)
      CLHEP::HepRandom::setTheSeed(time(0) + 42);
    
    // Only include visualization capabilities if Geant4 has been
    // built with the G4VIS_USE flag
    G4VisManager *visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // Color particle trajectories by particle type
    G4TrajectoryDrawByParticleID *colorModel = new G4TrajectoryDrawByParticleID;
    colorModel->Set("neutron", "cyan");
    colorModel->Set("gamma", "green");
    colorModel->Set("e-", "red");
    colorModel->Set("e+", "blue");
    colorModel->Set("proton", "yellow");
    colorModel->SetDefault("gray");
    visManager->RegisterModel(colorModel);
    visManager->SelectTrajectoryModel(colorModel->Name());
    
    // At present, two options are provided for visualization and
    // control of ZK with sequential architecture.  First, the hotness
    // of G4UI graphical user interface using Qt; second standard
    // command-line-plus-OpenGL
      if(visQt){
      // Build the G4UI GUI and run the Qt visualization macro
         G4UIExecutive *UIexecutive = new G4UIExecutive(argc, argv);
         UImanager->ApplyCommand("/control/execute {run}/ZK.Qt.vis");
         UImanager->ApplyCommand("/vis/scene/add/trajectories");
         UImanager->ApplyCommand("/vis/scene/add/hits");
         UIexecutive->SessionStart();
         delete UIexecutive;
      }
      else if(manualp>0)
      {
      // Changing to false will use same seed for all nodes
      if(true)
      CLHEP::HepRandom::setTheSeed(time(0) + (manualp*100));

      // ZK with MPI is setup to run all commands in batch mode from a
      // script. The default script is "runtime/ZK.mpi.mac" although the
      // user may specify a different macro as the 2nd cmd line arg.
      G4String macroCmd = "/control/execute " + MPIMacroName;
      UImanager->ApplyCommand(macroCmd);
      }
      else{
         if(visualization){
	         // Run the OpenGL visualization macro
	         UImanager->ApplyCommand("/control/execute {run}/ZK.OGLIX.vis");
	         UImanager->ApplyCommand("/vis/scene/add/trajectories");
	         UImanager->ApplyCommand("/vis/scene/add/hits");
         }
    
         // Create a decent 'tcsh'-like prompt for tab completion, command
         // history, etc.  Also, style points for cooler prompt
         G4String prompt = "ZK >> ";
         G4int maxHist = 200;
         G4UIsession* session = new G4UIterminal(new G4UItcsh(prompt, maxHist));


         // As Gallagher said: "Styyyyyyyyyyyle!"
         G4cout << "\n\n"
	        << "\t\t    ZZZZZZZ    K    K             \n"
	        << "\t\t          Z    K   K              \n"
	        << "\t\t         Z     K  K               \n"
	        << "\t\t        Z      KKK                \n"
    	    << "\t\t     ZZZZZ     K  K               \n"
	        << "\t\t      Z        K   K              \n"
	        << "\t\t     Z  ERO    K    K NOWLEDGE    \n"
	        << "\t\t    Z          K    K             \n"
	        << "\t\t    ZZZZZZ     K    K             \n"
	        << "\n\n      *******    WELCOME TO THE ZKBrem SIMULATION    *******\n\n\n"
	        << G4endl;
         
         session->SessionStart();
         delete session;
    
         if(visualization)
            delete visManager;
      }
      
  }
  
  /////////////////////////////////////
  // Current ZK instance is parallel //
  /////////////////////////////////////
  
#ifdef ZK_MPI_ENABLED  
  else{
    // Internally assign niceness of 19 to parallel builds of ZK
    G4int priority = 19;
    setpriority(PRIO_PROCESS,getpid(),priority);

    // The MPIManager requires explicit command line inputs.  Since ZK
    // needs to be configured by the user via the argv array, create
    // special argc/argv inputs for MPI parallelization
    //const G4int argcMPI = 2;
    //char *argvMPI[argcMPI];
    //argvMPI[0] = argv[0]; // binary name
    //argvMPI[1] = (char *)"/tmp/ZKSlave"; // slave file base name
    
    // Create the Open MPI manager
    //MPIManager *theMPIManager= new MPIManager(argcMPI,argvMPI);

    G4int MPI_Rank = theMPIManager->GetRank();

    // Changing to false will use same seed for all nodes
    if(true)
      CLHEP::HepRandom::setTheSeed(time(0) + 7*MPI_Rank);
    
    // ZK with MPI is setup to run all commands in batch mode from a
    // script. The default script is "runtime/ZK.mpi.mac" although the
    // user may specify a different macro as the 2nd cmd line arg.
    G4String macroCmd = "/control/execute " + MPIMacroName;
    UImanager->ApplyCommand(macroCmd);
    
    delete theMPIManager;
  }
#endif

  // General garbage collection

  delete theRSManager;

  //delete theMetadata;
  
  delete theRunManager;

  G4cout << "\n" << G4endl;
  
  return 0;
}

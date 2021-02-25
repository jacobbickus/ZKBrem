/////////////////////////////////////////////////////////////////////////////////
//
// name: MPIManager.hh
// date: 17 Sep 14
// auth: Zach Hartwig
// mail: hartwig@psfc.mit.edu
//
// desc: The MPIManager class handles most of the parallelization
//       functionality for custom builds of Geant4 simulations using
//       an instance of MPI (Open MPI, MPICH2, etc). It is constructed
//       as a meyer's singleton to provide a single instance of the
//       class that is available globally throughout the simulation
//       codebase (identical to the G4 manager concept). It is
//       responsible for such important functions as node
//       identification and control, event distribution, and data
//       reduction (i.e. aggregation of data contained on the slave
//       nodes to the master node for persistent output). In addition,
//       it provides a replacement method for launching particles for
//       parallel binaries that provides control of how the events are
//       distributed across the nodes as well as ability launch far
//       more events per run than the standard G4 /run/beamOn method.
//
//       Note that much of the MPIManager class code is protected by
//       compiler macros to enable non-MPI builds of the G4
//       simluation. Such macros can be used throughout the simulation
//       code to prevent compilation errors when MPI is not present on
//       the system.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef MPIManager_hh
#define MPIManager_hh 1

#include "globals.hh"

#include <fstream>

class MPIManagerMessenger;

class MPIManager
{
public:
  MPIManager(int, char **);
  ~MPIManager();

  // Static method to obtain the singleton MPI manager
  static MPIManager *GetInstance();

  // Method to run the beam used by MPImessenger
  void BeamOn(G4double, G4bool);

  // Methods to create and distribute seeds to nodes
  void CreateSeeds();
  void DistributeSeeds();

  // Method to obtain number of procs
  G4int GetSize() {return size;}

  // Method to obtain the node rank
  G4int GetRank() {return rank;}

  // Method for obtaining/operatong on the events to be processed
  G4double GetTotalEvents() {return totalEvents;}
  G4int GetMasterEvents() {return masterEvents;}
  G4int GetSlaveEvents() {return slaveEvents;}
  void ThrowEventError();

  // Method to determine whether events are distributed
  G4bool GetDistributeEvents() {return distributeEvents;}

  // Method to force an MPI barrier for all nodes
  void ForceBarrier(G4String);

  // Method to sum doubles on all nodes to a single value on the master
  G4double SumDoublesToMaster(G4double);

  // Method to sum ints on all nodes to a single value on the master
  G4int SumIntsToMaster(G4int);

private:
  // The singleton object
  static MPIManager *theMPImanager;

  // The messenger interface for run time commands
  MPIManagerMessenger *theMPImessenger;

  // size == number of procs; rank == process id
  G4int size, rank;
  G4bool isMaster, isSlave;
  enum {RANK_MASTER, RANK_SLAVE};
  std::ofstream slaveOut;

  G4double totalEvents;
  G4long masterEvents, slaveEvents;
  G4bool distributeEvents;

  std::vector<G4long> seedPacket;
};

#endif

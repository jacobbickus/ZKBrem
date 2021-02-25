#ifndef geometryConstruction_hh
#define geometryConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "sensitiveDetector.hh"
#include "G4SDManager.hh"
#include <vector>
using std::vector;
//#include "particleSplittingOperator.hh"

class geometryConstruction : public G4VUserDetectorConstruction
{
  
public:
  geometryConstruction(G4int ilysoconf = 0 , G4int iplateconf = 0, G4bool incff = false, G4bool incfilt = true);
  ~geometryConstruction();
  
  G4VPhysicalVolume *Construct();

  void BuildTarget(G4String name, vector<G4Material*> v_mats, vector<double> v_thickness, double width, double length, G4String opt = "D");

  G4LogicalVolume * BuildHPGe(unsigned int ind, G4ThreeVector pos, G4RotationMatrix * rot, double gA, double gB, double gC, double gD, double gE, double gF, double gG, double gH, double gI, double gJ, double gK, double gL, double gM, double gN);

private:

  G4int plateconf, lysoconf;
  G4bool incfilters, incfalsedets;
  G4double collimator_lead_z;

  G4Isotope *U235;
  G4Isotope *U238;
  G4Material *Uenriched95;
  G4Material *HMX;
  G4Material *WgPu_simple;

  G4Material *Vacuum;
  G4Material *Air;
  G4Material *Cu;
  G4Material *Au;
  G4Material *Al;
  G4Material *W;

  G4Material *Concrete;

  G4Material *Ge;
  G4Material *Pb;
  G4Material *AltPb;

  G4Material *B;
  G4Material *N;
  G4Material *C;
  G4Material *H;
  G4Material *O;
  G4Material *Si;

  G4Material *Cr;
  G4Material *Mn;
  G4Material *Fe;
  G4Material *Ni;

  G4Material *La;
  G4Material *Br;

  G4Material *U;
  G4Material *DU_mat;
  G4Material *Al_mat;

  G4Material *Lu;
  G4Material *Y;
  G4Material *Ce;

  G4Material *Al_alloy;
  G4Material *Mg;
  G4Material *Zn;
  G4Material *Ti;
  G4Material *Sn;

  G4double a;
  G4int iz, n;
  G4String name;
  G4String symbol;
  G4int n_components;
  G4double density, abundance, fraction_mass;

  G4Material *StainlessSteel;
  G4Material *Mylar;
  G4Material *LYSO_scint;
  G4Material *LYSO;
  G4Material *LaBr3_scint;
  G4Material *LaBr3;
  G4Material *BoratedPoly;

  G4VisAttributes *UVisAtt;
  G4VisAttributes *BPVisAtt;
  G4VisAttributes *WVisAtt;

  G4Box *world_S;
  G4LogicalVolume *world_L;
  G4VPhysicalVolume *world_P;

  G4Box *hall_S;
  G4LogicalVolume *hall_L;
  G4VPhysicalVolume *hall_P;

  G4Box *target_S;
  G4LogicalVolume *target_L;
  G4VPhysicalVolume *target_P;

  G4Tubs *hpge_S;
  G4LogicalVolume *hpge_L;
  G4VPhysicalVolume *hpge_P;

  G4Box *Cu_S;
  G4LogicalVolume *Cu_L;
  G4VPhysicalVolume *Cu_P;

  // G4Box *Au_S;
  G4LogicalVolume *Au_L;
  G4VPhysicalVolume *Au_P;

  G4Tubs *ss_S;
  G4LogicalVolume *ss_L;
  G4VPhysicalVolume *ss_P;

  G4Cons *cone_S;
  G4LogicalVolume *cone_L;
  G4VPhysicalVolume *cone_P;
  
  G4Box * lyso_S;
  G4LogicalVolume * lyso_L;
  G4VPhysicalVolume * lyso_P;
  
  G4VPhysicalVolume * vac_P;
  G4LogicalVolume * screw_L[4];
  G4VPhysicalVolume * screw_P[4];
  G4LogicalVolume * junk_L;
  G4LogicalVolume * beampipe_L;

  G4LogicalVolume * HPGe_crystal_L1;
  G4LogicalVolume * HPGe_crystal_L2;
  G4LogicalVolume * HPGe_crystal_L3;

  sensitiveDetector *hpgeSD_1;
  sensitiveDetector *hpgeSD_2;
  sensitiveDetector *hpgeSD_3;
  sensitiveDetector *NearSD;
  sensitiveDetector *FarSD;
  sensitiveDetector *lysoSD;

  G4double U_plate_z;

  G4bool useSplitting;
  // particleSplittingOperator* biasingOperator;
};

#endif

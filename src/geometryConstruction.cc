#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "sensitiveDetector.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trap.hh"
#include "G4GeometryManager.hh"
#include "G4PVReplica.hh"
#include "G4PSPassageCellCurrent3D.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSCellCharge.hh"

#include "geometryConstruction.hh"

//#include "G4GDMLParser.hh"


geometryConstruction::geometryConstruction(G4int ilysoconf, G4int iplateconf, G4bool incff, G4bool incfilt)
{
  lysoconf = ilysoconf;
  plateconf = iplateconf;
  incfilters = incfilt;
  incfalsedets = incff;
  collimator_lead_z = 0;

  G4NistManager *NISTMgr = G4NistManager::Instance();

  Vacuum = NISTMgr->FindOrBuildMaterial("G4_Galactic");
  Air = NISTMgr->FindOrBuildMaterial("G4_AIR");
  Cu = NISTMgr->FindOrBuildMaterial("G4_Cu");
  Au = NISTMgr->FindOrBuildMaterial("G4_Au");
  Al = NISTMgr->FindOrBuildMaterial("G4_Al");
  W  = NISTMgr->FindOrBuildMaterial("G4_W");

  Concrete = NISTMgr->FindOrBuildMaterial("G4_CONCRETE");

  Ge = NISTMgr->FindOrBuildMaterial("G4_Ge");
  Pb = NISTMgr->FindOrBuildMaterial("G4_Pb");

  C  = NISTMgr->FindOrBuildMaterial("G4_C");
  O  = NISTMgr->FindOrBuildMaterial("G4_O");
  H  = NISTMgr->FindOrBuildMaterial("G4_H");
  N  = NISTMgr->FindOrBuildMaterial("G4_N");
  B  = NISTMgr->FindOrBuildMaterial("G4_B");
  Si = NISTMgr->FindOrBuildMaterial("G4_Si");
  Cr = NISTMgr->FindOrBuildMaterial("G4_Cr");
  Mn = NISTMgr->FindOrBuildMaterial("G4_Mn");
  Fe = NISTMgr->FindOrBuildMaterial("G4_Fe");
  Ni = NISTMgr->FindOrBuildMaterial("G4_Ni");

  Mg = NISTMgr->FindOrBuildMaterial("G4_Mg");
  Si = NISTMgr->FindOrBuildMaterial("G4_Si");
  Zn = NISTMgr->FindOrBuildMaterial("G4_Zn");
  Ti = NISTMgr->FindOrBuildMaterial("G4_Ti");

  U = NISTMgr->FindOrBuildMaterial("G4_U");

  Lu = NISTMgr->FindOrBuildMaterial("G4_Lu");
  Y = NISTMgr->FindOrBuildMaterial("G4_Y");
  Ce = NISTMgr->FindOrBuildMaterial("G4_Ce");

  La = NISTMgr->FindOrBuildMaterial("G4_La");
  Br = NISTMgr->FindOrBuildMaterial("G4_Br");

  AltPb = new G4Material(name="AltPb", density=11.34*g/cm3, n_components = 2);
  AltPb->AddMaterial(Pb, fraction_mass=99.999*perCent);
  AltPb->AddMaterial(W, fraction_mass=0.001*perCent);

  LYSO_scint = new G4Material(name="LYSO_scintillator", density=7.4*g/cm3, n_components = 4);
  LYSO_scint->AddMaterial(Lu, fraction_mass=74.0*perCent);
  LYSO_scint->AddMaterial(Y, fraction_mass=2.0*perCent);
  LYSO_scint->AddMaterial(Si, fraction_mass=6.0*perCent);
  LYSO_scint->AddMaterial(O, fraction_mass=18.0*perCent);

  LYSO = new G4Material(name="LYSO_crystal", density=7.4*g/cm3, n_components = 2);
  LYSO->AddMaterial(LYSO_scint, fraction_mass=99.5*perCent);
  LYSO->AddMaterial(Ce, fraction_mass=0.5*perCent);

  LaBr3_scint = new G4Material(name="LaBr3_scintillator", density=5.08*g/cm3, n_components = 2);
  LaBr3_scint->AddMaterial(La, fraction_mass=36.7*perCent);
  LaBr3_scint->AddMaterial(Br, fraction_mass=63.3*perCent);

  LaBr3 = new G4Material(name="LaBr3_Ce", density=5.08*g/cm3, n_components = 2);
  LaBr3->AddMaterial(LaBr3_scint, fraction_mass=99.5*perCent);
  LaBr3->AddMaterial(Ce, fraction_mass=0.5*perCent);

  StainlessSteel = new G4Material(name="StainlessSteel", density= 8.06*g/cm3, n_components=6);
  StainlessSteel->AddMaterial(C, fraction_mass=0.001);
  StainlessSteel->AddMaterial(Si, fraction_mass=0.007);
  StainlessSteel->AddMaterial(Cr, fraction_mass=0.18);
  StainlessSteel->AddMaterial(Mn, fraction_mass=0.01);
  StainlessSteel->AddMaterial(Fe, fraction_mass=0.712);
  StainlessSteel->AddMaterial(Ni, fraction_mass=0.09);

  Mylar = new G4Material(name="Mylar", density= 1.4*g/cm3, n_components=3);
  Mylar->AddMaterial(H, fraction_mass=0.041959);
  Mylar->AddMaterial(C, fraction_mass=0.625017);
  Mylar->AddMaterial(O, fraction_mass=0.333025);

  // Borated Polyethylene
  BoratedPoly = new G4Material(name = "BoratedPoly", density = 1.000*g/cm3, n_components = 3);
  BoratedPoly->AddMaterial(H, fraction_mass=0.125355);
  BoratedPoly->AddMaterial(B, fraction_mass=0.100000);
  BoratedPoly->AddMaterial(C, fraction_mass=0.774645);

  // simplified DU that can actually undergo NRF
  U235 = new G4Isotope(name = "U235", iz = 92, n = 235, a = 235.044*g/mole);
  U238 = new G4Isotope(name = "U238", iz = 92, n = 238, a = 238.051*g/mole);
  G4Element *DU = new G4Element(name = "Depleted Uranium", symbol = "DU", n_components = 2);
  DU->AddIsotope(U235, abundance = 0.003);
  DU->AddIsotope(U238, abundance = 0.997);

  DU_mat = new G4Material(name = "Depleted Uranium", density = 18.95*g/cm3, n_components = 1);
  DU_mat->AddElement(DU, fraction_mass = 100.0*perCent);

  // Aluminum that can undergo NRF
  G4Isotope *Al27  = new G4Isotope(name = "Al27", iz = 13, n = 27, a = 26.982*g/mole);
  G4Element *NatAl = new G4Element(name = "Natural Aluminum", symbol = "NatAl", n_components = 1);
  NatAl->AddIsotope(Al27, abundance = 100*perCent);
  Al_mat = new G4Material(name = "Natural Aluminum", density = 2.6989*g/cm3, n_components = 1);
  Al_mat->AddElement(NatAl, fraction_mass = 100.0*perCent);

  Al_alloy = new G4Material(name = "Aluminum Alloy", density = 2.6989*g/cm3, n_components = 7);
  Al_alloy->AddMaterial(Al_mat, fraction_mass = 0.9729);
  Al_alloy->AddMaterial(Mg,     fraction_mass = 0.0136);
  Al_alloy->AddMaterial(Si,     fraction_mass = 0.0090);
  Al_alloy->AddMaterial(Cu,     fraction_mass = 0.0023);
  Al_alloy->AddMaterial(Fe,     fraction_mass = 0.0009);
  Al_alloy->AddMaterial(Zn,     fraction_mass = 0.0010);
  Al_alloy->AddMaterial(Ti,     fraction_mass = 0.0003);

  // HMX (high explosive)
  HMX = new G4Material(name = "HMX", density = 1.890*g/cm3, n_components = 4);
  HMX->AddMaterial(H, fraction_mass=0.027227);
  HMX->AddMaterial(C, fraction_mass=0.162222);
  HMX->AddMaterial(N, fraction_mass=0.378361);
  HMX->AddMaterial(O, fraction_mass=0.432190);

  G4Isotope *Pu239 = new G4Isotope(name = "Pu239", iz = 94, n = 239, a = 239.052*g/mole);
  G4Isotope *Pu240 = new G4Isotope(name = "Pu240", iz = 94, n = 240, a = 240.054*g/mole);

  // Simple uranium-238
  G4Element *PureU238 = new G4Element(name = "PureU238", symbol = "PureU238", n_components = 1);
  PureU238->AddIsotope(U238, abundance = 1.00);
  G4Material *U238_simple = new G4Material(name = "Simple U-238", density = 19.052*g/cm3, n_components = 1);
  U238_simple->AddElement(PureU238, fraction_mass = 100.0*perCent);

  // Simple uranium-235
  G4Element *PureU235 = new G4Element(name = "PureU235", symbol = "PureU235", n_components = 1);
  PureU235->AddIsotope(U235, abundance = 1.00);
  G4Material *U235_simple = new G4Material(name = "Simple U-235", density = 18.811*g/cm3, n_components = 1);
  U235_simple->AddElement(PureU235, fraction_mass = 100.0*perCent);

  Uenriched95 = new G4Material(name = "95pcw enriched uranium", density = 18.7*g/cm3, n_components = 2);
  Uenriched95->AddElement(PureU238, fraction_mass =  5.0*perCent);
  Uenriched95->AddElement(PureU235, fraction_mass = 95.0*perCent);

  // Simple plutonium-239
  G4Element *PurePu239 = new G4Element(name = "PurePu239", symbol = "PurePu239", n_components = 1);
  PurePu239->AddIsotope(Pu239, abundance = 1.00);
  G4Material *Pu239_simple = new G4Material(name = "Simple Pu-239", density = 19.41*g/cm3 , n_components = 1);
  Pu239_simple->AddElement(PurePu239, fraction_mass = 100.0*perCent);

  // Simple plutonium-240
  G4Element *PurePu240 = new G4Element(name = "PurePu240", symbol = "PurePu240", n_components = 1);
  PurePu240->AddIsotope(Pu240, abundance = 1.00);
  G4Material *Pu240_simple = new G4Material(name = "Simple Pu-240", density = 19.496*g/cm3 , n_components = 1);
  Pu240_simple->AddElement(PurePu240, fraction_mass = 100.0*perCent);

  WgPu_simple = new G4Material(name = "Simple WgPu", density = 19.84*g/cm3, n_components = 2);
  WgPu_simple->AddElement(PurePu239, fraction_mass = 0.94);
  WgPu_simple->AddElement(PurePu240, fraction_mass = 0.06);


  // Sensitive detector(s)
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();

  hpgeSD_1 = new sensitiveDetector("hpgeSD_1");
  SDMan->AddNewDetector(hpgeSD_1);

  hpgeSD_2 = new sensitiveDetector("hpgeSD_2");
  SDMan->AddNewDetector(hpgeSD_2);

  hpgeSD_3 = new sensitiveDetector("hpgeSD_3");
  SDMan->AddNewDetector(hpgeSD_3);

  lysoSD = new sensitiveDetector("lysoSD");
  SDMan->AddNewDetector(lysoSD);

  //NearSD = new sensitiveDetector("NearSD");
  //SDMan->AddNewDetector(NearSD);

//  FarSD = new sensitiveDetector("FarSD");
//  SDMan->AddNewDetector(FarSD);

  // particleSplittingOperator to attach to logical volumes
  //useSplitting = true;
  //if (useSplitting) biasingOperator = new particleSplittingOperator();
}


geometryConstruction::~geometryConstruction()
{;}


G4VPhysicalVolume * geometryConstruction::Construct()
{
  // ------------------------------------------------------------------------------------
  // The world volume
  //

  const G4double inch = 25.4; // one inch in (Geant4 default) mm, because lead bricks are 2x4x8"

  // distances from center of radiator to floor, ceiling
  G4double rad_to_floor = 110.5*cm;
  G4double rad_to_ceil  = 134.0*cm;

  // distances from back, front of radiator to back, front walls: need to account for radiator thickness
  G4double rad_to_back  = 330.0*cm;
  G4double rad_to_front = 445.0*cm;
  G4double rad_total = 20.1*mm;

  // distances from center of radiator to left, right walls (facing downbeam)
  G4double rad_to_left  = 382.0*cm;
  G4double rad_to_right = 259.0*cm;

  // assumed wall/ceiling/floor thickness of 30 cm
  G4double wall_thickness = 30.0*cm;

  G4double hallX = rad_to_left  + rad_to_right;
  G4double hallY = rad_to_floor + rad_to_ceil;
  G4double hallZ = rad_to_back  + rad_to_front + rad_total; // beam direction

  G4double worldX = hallX + 2.0*wall_thickness;
  G4double worldY = hallY + 2.0*wall_thickness;
  G4double worldZ = hallZ + 2.0*wall_thickness;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldZ);

  /*world_S = new G4Box("world_S", worldX/2, worldY/2, worldZ/2);

  world_L = new G4LogicalVolume(world_S,
        Concrete,
        "world_L");

  world_P = new G4PVPlacement(0,
            G4ThreeVector(),
            world_L,
            "TheWorld",
            0,
            false,
            0);*/

  G4VisAttributes *worldVisAtt = new G4VisAttributes();
  worldVisAtt->SetVisibility(false);

  hall_S = new G4Box("hall_S", hallX/2.0, hallY/2.0, hallZ/2.0);

  G4Material * worldmat = Air;
  if (!incfilters)
    worldmat = Vacuum;

  hall_L = new G4LogicalVolume(hall_S, worldmat, "hall_L");

  hall_P = new G4PVPlacement(0,G4ThreeVector(),hall_L,"TheHall",0/*world_L*/,false,0);

  hall_L->SetVisAttributes(worldVisAtt);

  // ------------------------------------------------------------------------------------
  // The bremsstrahlung target
  //
  bool oldgeo =  false; // Reimplementing the radiator geometry from scractch
                        // based on the measurements made at HVRL on 1-25-2017.
                        // Use "oldgeo" to get the old configuration.

  G4double CuX = 100.*mm; // Variables needed elsewhere
  G4double CuY = 100.*mm;
  G4double CuZ = 20.1*mm;

  if (oldgeo)
  {
    // Cu backing
    // all the radiator geometry is contained by the copper, so use it as mother volume

    // keep the radiator centered for now, since it's easier
    G4double radPosX = 0.0; //-worldX/2.0 + rad_to_right;
    G4double radPosY = 0.0; //-worldY/2.0 + rad_to_floor;
    G4double radPosZ = 0.0; //-worldZ/2.0 + rad_to_back;

    Cu_S = new G4Box("Cu_S", CuX/2.0, CuY/2.0, CuZ/2.0);

    Cu_L = new G4LogicalVolume(Cu_S,
           Cu,
           "Cu_L");

    Cu_P = new G4PVPlacement(0,
               G4ThreeVector(radPosX,radPosY,radPosZ),
               Cu_L,
               "CuPlate",
               hall_L,
               false,
               0);

    G4VisAttributes *CuVisAtt = new G4VisAttributes(G4Color(184./256,
                      115./256,
                      51./256,
                      1.0));
    CuVisAtt->SetForceSolid(false);
    Cu_L->SetVisAttributes(CuVisAtt);

     // Au disk  - 1 cm diameter, 102 um thickness

    G4double AuInR = 0.*mm;
    G4double AuOutR = 5.0*mm;
    G4double AuZ = 0.102*mm;
    G4double AuStartAng = 0.*degree;
    G4double AuSpanAng = 360.*degree;

    G4double AuPosZ = -(AuZ/2.0);

    G4Tubs * Au_S = new G4Tubs("Au_S", AuInR, AuOutR, AuZ/2.0, AuStartAng, AuSpanAng);

    Au_L = new G4LogicalVolume(Au_S,
           Au,
           "Au_L");

    Au_P = new G4PVPlacement(0,
               G4ThreeVector(0., 0., AuPosZ),
               Au_L,
               "TheGoldLayer",
               Cu_L,
               false,
               0);

    G4VisAttributes *AuVisAtt = new G4VisAttributes(G4Color(243./256,
                      214./256,
                      26./256,
                      1.0));
    AuVisAtt->SetForceSolid(false);
    Au_L->SetVisAttributes(AuVisAtt);

    // steel plate in front of gold layer - 4 cm diameter, 1 cm thickness

    G4double ssInR = 0.*mm;
    G4double ssOutR = 30.*mm;
    G4double ssZ = (CuZ/2-AuZ)*mm;
    G4double ssStartAng = 0.*degree;
    G4double ssSpanAng = 360.*degree;

    G4double ssPosZ = -(AuZ+ssZ/2.0)*mm;

    ss_S = new G4Tubs("ss_S", ssInR, ssOutR, ssZ/2.0, ssStartAng, ssSpanAng);

    ss_L = new G4LogicalVolume(ss_S,
           StainlessSteel,
           "ss_L");

    ss_P = new G4PVPlacement(0,
               G4ThreeVector(0., 0., ssPosZ),
               ss_L,
               "TheStainlessSteelPlate",
               Cu_L,
               false,
               0);

    G4VisAttributes *ssVisAtt = new G4VisAttributes(G4Color(115./256,
                      115./256,
                      115./256,
                      0.5));
    ssVisAtt->SetForceSolid(false);
    ss_L->SetVisAttributes(ssVisAtt);

    // truncated cone hole inside stainless steel plate

    G4double  coneRmin1 = 0.*mm;
    G4double  coneRmax1 = 10.*mm;
    G4double  coneRmin2 = 0.*mm;
    G4double  coneRmax2 = 5.0*mm;
    G4double  coneZ = (CuZ/2-AuZ)*mm;
    G4double  coneSPhi = 0.*degree;
    G4double  coneDPhi = 360.*degree;

    //G4double conePosZ = -(AuZ+coneZ/2.0)*mm;

    cone_S = new G4Cons("cone_S", coneRmin1, coneRmax1, coneRmin2, coneRmax2, coneZ/2.0, coneSPhi, coneDPhi);

    cone_L = new G4LogicalVolume(cone_S,
           Vacuum,
           "cone_L");

    cone_P = new G4PVPlacement(0,
               G4ThreeVector(0., 0., 0.),
               cone_L,
               "TheConicalHole",
               ss_L,
               false,
               0);

    cone_L->SetVisAttributes(ssVisAtt);
  }
  else // Redone radiator geometry
  {

    // Main Cu block with recesses of vacuum, implemented via booleans since the
    // foil mount disc extrudes from the recess, making containing solids a pain

    CuX = 120.65*mm;
    CuY = 120.65*mm;
    CuZ = 25.4*mm;

    // Will need to adjust away from zero pos as construction is compelte
    G4double radPosX = 0.0; //-worldX/2.0 + rad_to_right;
    G4double radPosY = 0.0; //-worldY/2.0 + rad_to_floor;
    G4double radPosZ = 0.0; //-worldZ/2.0 + rad_to_back;

    Cu_S = new G4Box("Cu_S", CuX/2.0, CuY/2.0, CuZ/2.0);

    // Wide/shallow recess
    G4double drec1 = 52.3*mm;
    G4double hrec1 = 2.48*mm;

    G4Tubs * rec1_S = new G4Tubs("Recess1_S", 0.0, drec1/2.0, hrec1/2.0, 0.0, 2.0*M_PI);
    G4SubtractionSolid * Cu1_S = new G4SubtractionSolid("CuCut1",Cu_S,rec1_S,0,G4ThreeVector(0,0,-CuZ/2+hrec1/2.0));


    // Deep, post-foil recess
    G4double drec2 = 11.05*mm;
    G4double hrec2 = 12.7*mm;

    G4Tubs * rec2_S = new G4Tubs("Recess2_S", 0.0, drec2/2.0, hrec2/2.0, 0.0, 2.0*M_PI);
    G4SubtractionSolid * Cu2_S = new G4SubtractionSolid("CuCut2",Cu1_S,rec2_S,0,G4ThreeVector(0,0,-CuZ/2+hrec1+hrec2/2.0));

    // Pressure balance channel
    G4Tubs * rec3_S = new G4Tubs("Recess3_S", 0.0, 2*mm/2.0, (16.4*mm*2+drec2)/2.0, 0.0, 2.0*M_PI);
    G4RotationMatrix * chanrot = new G4RotationMatrix(); chanrot->rotateY(M_PI/2.0);
    G4SubtractionSolid * Cu3_S = new G4SubtractionSolid("CuCut3",Cu2_S,rec3_S,chanrot,G4ThreeVector(0,0,-CuZ/2+hrec1));

    Cu_L = new G4LogicalVolume(Cu3_S,Cu,"Cu_L");

    // Use these dimensions to create a vacuum volume to contain the entirety of
    // the radiator apparatus and represent the inside of the beam pipe
    G4Tubs * beamvac_S = new G4Tubs("beamVacuum_S",0.0,drec1/2.0,20*cm,0,2.0*M_PI);
    G4UnionSolid * allvac_S = new G4UnionSolid("vacuum_S",Cu_S,beamvac_S,0,G4ThreeVector(0,0,-CuZ/2-20*cm));

    G4LogicalVolume * allvac_L =  new G4LogicalVolume(allvac_S,Vacuum,"allvac_L");

    vac_P = new G4PVPlacement(0,
           G4ThreeVector(radPosX,radPosY,radPosZ),
           allvac_L,
           "vacuumline",
           hall_L,
           false,
           0);

    Cu_P = new G4PVPlacement(0,
               G4ThreeVector(0,0,0),
               Cu_L,
               "CuPlate",
               allvac_L,
               false,
               0);

    // Gold foil
    G4double AuThickness = 0.126*mm;
    G4Box * Au_S = new G4Box("Au_S", 3.0/4.0*inch/2.0, 3.0/4.0*inch/2.0, AuThickness/2.0);

    Au_L = new G4LogicalVolume(Au_S,Au,"Au_L");

    double AuPosZ = -CuZ/2+hrec1-AuThickness/2;

    Au_P = new G4PVPlacement(0,
               G4ThreeVector(0., 0., AuPosZ),
               Au_L,
               "TheGoldLayer",
               allvac_L,
               false,
               0);

    // Deposited metal vapor on Au surface
    double junkT = 0.3*mm;
    G4Tubs * junk_S = new G4Tubs("junk_S", 0.0, 9.5*mm/2.0, junkT/2.0, 0.0, 2.0*M_PI);
    junk_L = new G4LogicalVolume(junk_S,StainlessSteel,"junk_L");
//    G4VPhysicalVolume * junk_P = new G4PVPlacement(0,
//                                   G4ThreeVector(0., 0., AuPosZ-AuThickness/2-junkT/2),
//                                   junk_L,
//                                   "Junk",
//                                   allvac_L,
//                                   false,
//                                   0);


    // Conical-hole disk foil holder thing (will still be called "ss_L", etc.
    // for the charge counter, even though we found out it's actually made of
    // copper)
    double diskh = 5.02*mm;
    G4Tubs * disk_S = new G4Tubs("disk_S", 0.0, 1.375*inch/2.0, diskh*mm/2.0, 0.0, 2.0*M_PI);
    G4Tubs * dbore1_S = new G4Tubs("dbore1_S",0.0, 10*mm/2.0, 10*mm, 0.0, 2.0*M_PI);
    G4SubtractionSolid * disk1_S = new G4SubtractionSolid("disk1_S",disk_S,dbore1_S,0,G4ThreeVector(0,0,0));
    cone_S = new G4Cons("cone_S", 0.0, 18.0*mm/2.0, 0.0, 10.0*mm/2.0, (diskh-1.0*mm)/2.0, 0.0, 2.0*M_PI);
    G4SubtractionSolid * disk2_S = new G4SubtractionSolid("disk2_S",disk1_S,cone_S,0,G4ThreeVector(0,0,-diskh/2.0+(diskh-1.0*mm)/2.0));

    ss_L = new G4LogicalVolume(disk2_S,Cu,"ss_L");
    ss_P = new G4PVPlacement(0,
           G4ThreeVector(0., 0., AuPosZ-AuThickness/2-diskh/2),
           ss_L,
           "DiskMount",
           allvac_L,
           false,
           0);

    // Disk mount screws
    G4Tubs * screw_S = new G4Tubs("screw_S", 0.0, 4.0*mm/2.0, 3.0*mm/2.0, 0.0, 2.0*M_PI);
    double screwx[4] = {9.23*mm,9.23*mm,-9.23*mm,-9.23*mm};
    double screwy[4] = {9.23*mm,-9.23*mm,9.23*mm,-9.23*mm};
    for (unsigned int j=0; j<4; j++)
    {
      char buff1[128];
      sprintf(buff1,"ss_L%d",j);
      screw_L[j] =  new G4LogicalVolume(screw_S,StainlessSteel,buff1);
      sprintf(buff1,"DiskScrew%d",j);
      screw_P[j] = new G4PVPlacement(0,
                        G4ThreeVector(screwx[j],screwy[j], AuPosZ-AuThickness/2-diskh-3.0*mm/2.0),
                        screw_L[j],
                        buff1,
                        allvac_L,
                        false,
                        0);
    }

    // Beam pipe section in electrical contact
    G4Tubs * beampipe_S = new G4Tubs("beampipe_out_S", drec1/2.0,drec1/2.0+3*mm, 2*inch/2, 0, 2.0*M_PI);
    beampipe_L = new G4LogicalVolume(beampipe_S, StainlessSteel, "beampipe_L");
    G4VPhysicalVolume * beampipe_P = new G4PVPlacement(0,
                               G4ThreeVector(0., 0., -CuZ/2.0-2.0*inch/2.0),
                               beampipe_L,
                               "ContactPipe",
                               hall_L,
                               false,
                               0);
  }

  // Charge accumulation scorers for the metal radiator components
  G4MultiFunctionalDetector * cuScore = new G4MultiFunctionalDetector("cuScore");
  G4MultiFunctionalDetector * auScore = new G4MultiFunctionalDetector("auScore");
  G4MultiFunctionalDetector * steelScore = new G4MultiFunctionalDetector("steelScore");
  G4SDManager::GetSDMpointer()->AddNewDetector(cuScore);
  G4SDManager::GetSDMpointer()->AddNewDetector(auScore);
  G4SDManager::GetSDMpointer()->AddNewDetector(steelScore);
  G4VPrimitiveScorer * cuChargeDep = new G4PSCellCharge("cuChargeDep");
  G4VPrimitiveScorer * auChargeDep = new G4PSCellCharge("auChargeDep");
  G4VPrimitiveScorer * steelChargeDep = new G4PSCellCharge("steelChargeDep");
  cuScore->RegisterPrimitive(cuChargeDep);
  auScore->RegisterPrimitive(auChargeDep);
  steelScore->RegisterPrimitive(steelChargeDep);

  G4MultiFunctionalDetector * junkScore = new G4MultiFunctionalDetector("junkScore");
  G4SDManager::GetSDMpointer()->AddNewDetector(junkScore);
  G4VPrimitiveScorer * junkChargeDep = new G4PSCellCharge("junkChargeDep");
  junkScore->RegisterPrimitive(junkChargeDep);

  if (junk_L)
    junk_L->SetSensitiveDetector(junkScore);

  if (!oldgeo)
  {
    for (unsigned int j=0; j<4; j++)
    {
        char buff1[128];
        sprintf(buff1,"screwScore%d",j);
        G4MultiFunctionalDetector * screwScore = new G4MultiFunctionalDetector(buff1);
        G4SDManager::GetSDMpointer()->AddNewDetector(screwScore);
        sprintf(buff1,"screwChargeDep%d",j);
        G4VPrimitiveScorer * screwChargeDep = new G4PSCellCharge(buff1);
        screwScore->RegisterPrimitive(screwChargeDep);
        screw_L[j]->SetSensitiveDetector(screwScore);
    }

    G4MultiFunctionalDetector * pipeScore = new G4MultiFunctionalDetector("pipeScore");
    G4SDManager::GetSDMpointer()->AddNewDetector(pipeScore);
    G4VPrimitiveScorer * pipeChargeDep = new G4PSCellCharge("pipeChargeDep");
    pipeScore->RegisterPrimitive(pipeChargeDep);
    beampipe_L->SetSensitiveDetector(pipeScore);

  }

  Cu_L->SetSensitiveDetector(cuScore);
  Au_L->SetSensitiveDetector(auScore);
  ss_L->SetSensitiveDetector(steelScore);

  const bool testbool = true;
  if (testbool){
  // ------------------------------------------------------------------------------------
  // Some useful lab dimensions
  //
  G4double lead_brick_2in = 2.0*inch;
  G4double lead_brick_4in = 4.0*inch;
  G4double lead_brick_8in = 8.0*inch;

  G4double lead_filter_thickness =  0.5*inch;
  G4double lead_filter_height    = 16.0*inch;
  G4double lead_filter_width     =  8.0*inch;

  G4double   rad_to_far_shield = 97.0*cm; // z distance from radiator front to center of far assembly
  G4double plane_to_far_shield =  2.0*cm; // y distance from beam plane to center height of far assembly
  G4double  beam_to_far_shield = 44.0*cm; // x distance from radiator center to center of far assembly at outer shield edge

  G4double rad_to_lyso = 63*inch + 1007.9; // Distance from radiator to LYSO crystal (not the same for all runs, distance was measured from plate)
  if (lysoconf == 0 ) rad_to_lyso = 96*inch + 1007.9;
  if (lysoconf == 3 ) rad_to_lyso = 143*inch - 2.43 + 11.0; // April 2017
  if (lysoconf == 4 ) rad_to_lyso = 143*inch - 2.43 + 0.75*inch;

  //
  // LYSO or LaBr3 Crystal
  //
  // LYSO crystal
  G4double lysoPosZ = rad_to_lyso;
  G4double lysoZ = 22.0*mm;
  if (lysoconf<4)
  {
    G4double lysoX = 4.0*mm;
    G4double lysoY = 4.0*mm;

    lyso_S = new G4Box("lyso_S", lysoX/2.0, lysoY/2.0, lysoZ/2.0);

    lyso_L = new G4LogicalVolume(lyso_S,
				   LYSO,
				   "lyso_L");

    lyso_P = new G4PVPlacement(0,
			         G4ThreeVector(0., 0., lysoPosZ),
			         lyso_L,
			         "TheLYSOcrystal",
			         hall_L,
			         false,
			         0);

    G4VisAttributes *lysoVisAtt = new G4VisAttributes(G4Color(184./256,
							      115./256,
							      51./256,
							      1.0));
    lysoVisAtt->SetForceSolid(true);

    lyso_L->SetSensitiveDetector(lysoSD);
  }
  else
  {
    lysoZ = 1.5*inch;

    G4Tubs * LaBr3_S = new G4Tubs("LaBr3_S", 0.0, 1.5*inch/2.0, 1.5*inch/2.0, 0, 2*M_PI);

    lyso_L = new G4LogicalVolume(LaBr3_S,
				   LaBr3,
				   "LaBr3_L");

    lyso_P = new G4PVPlacement(0,
			         G4ThreeVector(0., 0., lysoPosZ),
			         lyso_L,
			         "TheLaBr3crystal",
			         hall_L,
			         false,
			         0);

    G4VisAttributes *lysoVisAtt = new G4VisAttributes(G4Color(184./256,
							      115./256,
							      51./256,
							      1.0));
    lysoVisAtt->SetForceSolid(true);

    lyso_L->SetSensitiveDetector(lysoSD);
  }

  // LYSO shielding (3 possible configurations, and no shielding for other numbers)
  if (lysoconf==0 || lysoconf==1) // Single brick (earliest HVRL runs)
  {
    G4Box * ls0_S = new G4Box("ls0_S", lead_brick_4in/2.0, lead_brick_8in/2.0, lead_brick_2in/2.0);
    G4LogicalVolume * ls0_L = new G4LogicalVolume(ls0_S, Pb, "ls0_L");
    G4VPhysicalVolume * ls0_P = new G4PVPlacement(0,
			              G4ThreeVector(0., 0., lysoPosZ-lysoZ/2-lead_brick_2in/2.0-5*mm),
			              ls0_L,
			              "LYSOshield0",
			              hall_L,
			              false,
			              0);
  }
  if (lysoconf==1) // Add second brick (present during most U plate runs)
  {
    G4Box * ls1_S = new G4Box("ls1_S", lead_brick_4in/2.0, lead_brick_8in/2.0, lead_brick_2in/2.0);
    G4LogicalVolume * ls1_L = new G4LogicalVolume(ls1_S, Pb, "ls1_L");
    G4VPhysicalVolume * ls1_P = new G4PVPlacement(0,
			              G4ThreeVector(0., 0., lysoPosZ-lysoZ/2-lead_brick_2in/2.0-lead_brick_2in-5*mm),
			              ls1_L,
			              "LYSOshield1",
			              hall_L,
			              false,
			              0);
  }
  if (lysoconf==2) // "Hut" shielding configuration, present during most Al runs
  {
    G4Box * ls_S[3];
    G4LogicalVolume * ls_L[3];
    G4PVPlacement * ls_P[3];
    double xlsp[3] = {0 , +lead_brick_2in, -lead_brick_2in};
    double zlsp[3] = {2*cm+lead_brick_2in, 2*inch, 1.5*inch};
    for (unsigned int j=0; j<3; j++)
    {
      char buff[512];
      sprintf(buff,"ls_S%d",j);
      ls_S[j] = new G4Box(buff, lead_brick_2in/2.0, lead_brick_8in/2.0, lead_brick_4in/2.0);
      sprintf(buff,"ls_L%d",j);
      ls_L[j] = new G4LogicalVolume(ls_S[j], Pb, buff);
      sprintf(buff,"LYSOshield%d",j);
      ls_P[j] = new G4PVPlacement(0,
			                G4ThreeVector(xlsp[j], 0., lysoPosZ - lysoZ/2 - zlsp[j]),
			                ls_L[j],
			                buff,
			                hall_L,
			                false,
			                0);
    }
    G4Box * ls1_S = new G4Box("ls1_S", lead_brick_8in/2.0, lead_brick_2in/2.0, lead_brick_4in/2.0);
    G4LogicalVolume * ls1_L = new G4LogicalVolume(ls1_S, Pb, "ls1_L");
    G4VPhysicalVolume * ls1_P = new G4PVPlacement(0,
			              G4ThreeVector(0., lead_brick_4in+lead_brick_2in/2.0, lysoPosZ-lysoZ/2.0-2.0*inch),
			              ls1_L,
			              "LYSOshield1",
			              hall_L,
			              false,
			              0);
  }
  if (lysoconf==3 || lysoconf==4) // April 2017 configuration or September 2017 configuration with 2 more inches of Pb
  {
    G4Box * ls_S[8];
    G4LogicalVolume * ls_L[8];
    G4PVPlacement * ls_P[8];
    double ylsp[8] = {-5.0*inch, +5.0*inch, -5.0*inch, -3.0*inch, -1.0*inch, 1.0*inch, 3.0*inch, 5.0*inch};
    double fwall = rad_to_lyso-4.0*inch;
    double zlsp[8] = {rad_to_lyso,rad_to_lyso,fwall,fwall,fwall,fwall,fwall,fwall};
    double addthick = 0.0;

    for (unsigned int j=0; j<8; j++)
    {

      if (lysoconf==4 && j>1) // September 2017 extra shielding
        addthick = 1.0*inch;

      char buff[512];
      sprintf(buff,"ls_S%d",j);
      ls_S[j] = new G4Box(buff, lead_brick_8in/2.0, lead_brick_2in/2.0, lead_brick_4in/2.0+addthick);
      sprintf(buff,"ls_L%d",j);
      ls_L[j] = new G4LogicalVolume(ls_S[j], Pb, buff);
      sprintf(buff,"LYSOshield%d",j);
      ls_P[j] = new G4PVPlacement(0,
			                G4ThreeVector(0, ylsp[j], zlsp[j]-addthick),
			                ls_L[j],
			                buff,
			                hall_L,
			                false,
			                0);
    }
    G4Box * ls1_S = new G4Box("ls1_S", lead_brick_2in/2.0, lead_brick_8in/2.0, lead_brick_4in/2.0);
    G4LogicalVolume * ls1_L = new G4LogicalVolume(ls1_S, Pb, "ls1_L");
    G4VPhysicalVolume * ls1_P = new G4PVPlacement(0,
			              G4ThreeVector(-3.0*inch, 0.0, rad_to_lyso),
			              ls1_L,
			              "LYSOshield1",
			              hall_L,
			              false,
			              0);
    G4Box * ls2_S = new G4Box("ls2_S", lead_brick_2in/2.0, lead_brick_8in/2.0, lead_brick_4in/2.0);
    G4LogicalVolume * ls2_L = new G4LogicalVolume(ls2_S, Pb, "ls1_L");
    G4VPhysicalVolume * ls2_P = new G4PVPlacement(0,
			              G4ThreeVector(+3.0*inch, 0.0, rad_to_lyso),
			              ls1_L,
			              "LYSOshield2",
			              hall_L,
			              false,
			              0);
  }
  // ------------------------------------------------------------------------------------
  // HPGe detectors
  //

  // our GEM80 detector; specs from ORTEC tech support
  // is a solid cylinder here, not coax like in reality
  G4double dHPGe_det = 81.8*mm;       // crystal diameter
  // G4double rHPGe_det = dHPGe_det/2.0; // crystal radius
  G4double hHPGe_det = 69.4*mm;       // crystal height

  // More detailed detector implementation, using ORTEC spec sheets.  Variables
  // defined here correspond to the same ones as on the spec sheet so that you can
  // quickly switch which crystal you are using.

  // DET
  double gA = dHPGe_det; // Crystal diameter
  double gB = hHPGe_det; // Crystal length
  double gC = 8.7*mm; // Bore diameter
  double gD = 55.9*mm; // Bore depth
  double gE = 5*mm;  // Bore end radius
  double gF = 130*mm; // Cup length
  double gG = 4*mm; // Cup to aluminized Mylar
  double gH = 0.03*mm; // Aluminized Mylar thickness (each of Al and Mylar) (cup cap)
  double gI = 1.0*mm; // Aluminum housing cap thickness
  double gJ = 8.0*mm; // Crystal corner rounding radius
  double gK = 0.8*mm; // Cup side thickness (aluminum)
  double gL = 1.5*mm; // Aluminum housing side thickness
  double gM = 0.7*mm; // Outer dead layer (Ge/Li) thickness
  double gN = 0.3*um; // Inner dead layer (Ge/B) thickness

  // other two detectors coded up after the stainless steel housings are defined

  G4RotationMatrix *HPGe_rot_2 = new G4RotationMatrix();
  G4double rot_angle_2 = pi/2.0; // rotate to 90 deg
  HPGe_rot_2->rotateY(rot_angle_2);

  G4RotationMatrix *HPGe_rot_1 = new G4RotationMatrix();
  G4double relative_angle = 35*pi/180.0; // 35 deg wrt det2
  HPGe_rot_1->rotateY(rot_angle_2-relative_angle);

  G4RotationMatrix *nHPGe_rot_2 = new G4RotationMatrix();
  G4double nrot_angle_2 = -pi/2.0; // rotate to 90 deg
  nHPGe_rot_2->rotateY(nrot_angle_2);

  G4RotationMatrix *nHPGe_rot_1 = new G4RotationMatrix();
  G4double nrelative_angle = -35*pi/180.0; // 35 deg wrt det2
  nHPGe_rot_1->rotateY(nrot_angle_2-nrelative_angle);


  // nominal positions of the centers of the HPGe crystals if they were pushed right up against the lead shields
  // in reality, they are probably a bit further back in their housings, but these are still useful coordinates

  // centerline in the det stack is 2 cm above beam axis, then subtract half the housing height, then another inch
  G4double det2_h_above_base = 5.0*inch/2.0 + 1.0*inch;

  G4double det2x = beam_to_far_shield + lead_filter_thickness + hHPGe_det/2.0;
  G4double det2y = plane_to_far_shield - 2.5*inch; // det centered in 5" housing, but housing is 2 cm above beam axis
  G4double det2z = CuZ/2.0 + rad_to_far_shield;

  G4double det1x = det2x - hHPGe_det/2.0 - lead_filter_width/2.0 * sin(relative_angle) + hHPGe_det/2.0 * cos(relative_angle);
  G4double det1y = det2y;
  G4double det1z = det2z - lead_filter_width/2.0 - lead_filter_thickness - lead_brick_2in - lead_filter_width/2.0 * cos(relative_angle) - hHPGe_det/2.0 * sin(relative_angle);

  // instead of the nominal coordinates above, the detectors may be further back in their housings
  G4double det3x_actual = det1x;
  G4double det3y_actual = det1y;
  G4double det3z_actual = det1z;

  G4double det2x_actual = det1x;
  G4double det2y_actual = det1y + 2.5*lead_brick_2in;
  G4double det2z_actual = det1z;

  G4double det1x_actual = det1x;
  G4double det1y_actual = det1y;
  G4double det1z_actual = det1z;

  // How far back in the steel box each sits
  G4double det3_displacement = 2.6*inch-1.75*inch;
  if (det3_displacement != 0.0){
    det3x_actual = det3x_actual + det3_displacement * cos(relative_angle);
    det3y_actual = det3y_actual;
    det3z_actual = det3z_actual - det3_displacement * sin(relative_angle)-0.5*inch;;
  }

  G4double det2_displacement = 2.25*inch;
  if (det2_displacement != 0.0){
    det2x_actual = det2x_actual + det2_displacement * cos(relative_angle);
    det2y_actual = det2y_actual;
    det2z_actual = det2z_actual - det2_displacement * sin(relative_angle);;
  }

  G4double det1_displacement = 2.82*inch;
  if (det1_displacement != 0.0){
    det1x_actual = det1x_actual + det1_displacement * cos(relative_angle);
    det1y_actual = det1y_actual;
    det1z_actual = det1z_actual - det1_displacement * sin(relative_angle);

  }

  std::cout<<"\n\nDetector 1 box: "<<det1x_actual<<" "<<det1y_actual<<" "<<det1z_actual<<"\n\n";

  G4ThreeVector pos1(det1x_actual-6.3*sin(relative_angle)-5.0*cos(relative_angle) + 3.0*inch*cos(relative_angle), det1y_actual-6.3, det1z_actual- 6.3*cos(relative_angle) + 5.0*sin(relative_angle) -3.0*inch*sin(relative_angle));

  G4ThreeVector pos2(det2x_actual-6.3*sin(relative_angle)-5.0*cos(relative_angle) + 3.0*inch*cos(relative_angle), det2y_actual-6.3, det2z_actual- 6.3*cos(relative_angle) + 5.0*sin(relative_angle) -3.0*inch*sin(relative_angle));

  G4ThreeVector pos3(-(det3x_actual-6.3*sin(relative_angle)-5.0*cos(relative_angle) + 3.0*inch*cos(relative_angle)), det3y_actual-6.3, det3z_actual- 6.3*cos(relative_angle) + 5.0*sin(relative_angle) -3.0*inch*sin(relative_angle));

  G4LogicalVolume * hl1 = BuildHPGe(0, pos1, HPGe_rot_1, gA, gB, gC, gD, gE, gF, gG, gH, gI, gJ, gK, gL, 0.573, gN);
  HPGe_crystal_L1 = hl1;
  G4LogicalVolume * hl2 = BuildHPGe(1, pos2, HPGe_rot_1, 81.2, 85.8, 11.6, 72.5, gE, gF, gG, gH, 1.5, gJ, gK, gL, 0.112, gN);
  HPGe_crystal_L2 = hl2;
  G4LogicalVolume * hl3 = BuildHPGe(2, pos3, nHPGe_rot_1, 78.4, 82.3, 9.6, 69.0, gE, gF, gG, gH, gI, gJ, gK, gL, 4.04, gN);
  HPGe_crystal_L3 = hl3;


  G4VisAttributes *GeVisAtt = new G4VisAttributes(G4Color(87./256, 193./256, 84./256));
  GeVisAtt->SetVisibility(false);
  GeVisAtt->SetForceSolid(false);
  HPGe_crystal_L1->SetVisAttributes(GeVisAtt);
  HPGe_crystal_L2->SetVisAttributes(GeVisAtt);
  HPGe_crystal_L3->SetVisAttributes(GeVisAtt);

  HPGe_crystal_L1->SetSensitiveDetector(hpgeSD_1);
  HPGe_crystal_L2->SetSensitiveDetector(hpgeSD_2);
  HPGe_crystal_L3->SetSensitiveDetector(hpgeSD_3);

  // ------------------------------------------------------------------------------------
  // Lead shielding
  //

  // half-inch lead filters on the front faces of detectors
  // more useful for positioning than actually blocking bremsstrahlung

  G4Box *lead_filter_S = new G4Box("lead_filter_S", lead_filter_width/2.0, lead_filter_height/2.0, lead_filter_thickness/2.0);

  G4LogicalVolume *lead_filter_L = new G4LogicalVolume(lead_filter_S, Pb, "lead_filter_L");

  G4double filter2x = det2x - hHPGe_det/2.0 - lead_filter_thickness/2.0;
  G4double filter2y = det2y - det2_h_above_base + lead_filter_height/2.0;
  G4double filter2z = det2z;

  G4double filter1x = det1x - (lead_filter_thickness + hHPGe_det)/2.0 * cos(relative_angle);
  G4double filter1y = filter2y;
  G4double filter1z = det1z + (lead_filter_thickness + hHPGe_det)/2.0 * sin(relative_angle);

  // additional lead filter thickness on det1
  G4double extra_filter_thickness = 0.5*inch;
  G4double extra_filter_length = 8.0*inch;
  G4double extra_filter_height = 16.0*inch;
  G4Box *extra_filter_S = new G4Box("extra_filter_S", extra_filter_length/2.0, extra_filter_height/2.0, extra_filter_thickness/2.0);
  G4LogicalVolume *extra_filter_L = new G4LogicalVolume(extra_filter_S, AltPb, "extra_filter_L");

  G4double tlby = det2y - det2_h_above_base - lead_brick_2in/2.0;
  G4double extra_filter_x = filter1x - (extra_filter_thickness + lead_filter_thickness)/2.0 * cos(relative_angle);
  G4double extra_filter_y = tlby + lead_brick_2in/2.0 + extra_filter_height/2.0;
  G4double extra_filter_z = filter1z + (extra_filter_thickness + lead_filter_thickness)/2.0 * sin(relative_angle);

  if (incfilters)
  {
    G4PVPlacement *lead_filter_1_P = new G4PVPlacement(HPGe_rot_1,
                                                       G4ThreeVector(filter1x, filter1y, filter1z),
                                                       extra_filter_L,
                                                       "LeadFilter1",
                                                       hall_L,
                                                       false,
                                                       0);

    G4PVPlacement *lead_filter_2_P = new G4PVPlacement(HPGe_rot_2,
                                                       G4ThreeVector(filter2x, filter2y, filter2z),
                                                       lead_filter_L,
                                                       "LeadFilter2",
                                                       hall_L,
                                                       false,
                                                       0);

    G4PVPlacement *lead_filter_3_P = new G4PVPlacement(nHPGe_rot_1,
                                                       G4ThreeVector(-filter1x+1.43*inch, filter1y, filter1z),
                                                       extra_filter_L,
                                                       "LeadFilter3",
                                                       hall_L,
                                                       false,
                                                       0);

    G4PVPlacement *lead_filter_4_P = new G4PVPlacement(nHPGe_rot_2,
                                                       G4ThreeVector(-filter2x+1.43*inch, filter2y, filter2z),
                                                       lead_filter_L,
                                                       "LeadFilter4",
                                                       hall_L,
                                                       false,
                                                       0);
  }

  // lead shielding behind the radiator
  G4double lead_behind_offset = 2.0*cm; // FIXME best guess
  G4Box *lead_behind_radiator_box = new G4Box("lead_behind_radiator_box", 6*lead_brick_2in/2.0, lead_brick_8in/2.0, 1.5*lead_brick_4in/2.0);
  G4Box *lead_behind_radiator_cut = new G4Box("lead_behind_radiator_cut", 30*mm, 30*mm, 600);
  G4SubtractionSolid *lead_behind_radiator_S = new G4SubtractionSolid("lead_behind_radiator_S", lead_behind_radiator_box, lead_behind_radiator_cut);
  G4LogicalVolume *lead_behind_radiator_L = new G4LogicalVolume(lead_behind_radiator_S, Pb, "lead_behind_radiator_L");
  G4PVPlacement *lead_behind_radiator_P = new G4PVPlacement(0,
                                                G4ThreeVector(0,0, CuZ/2+3.0+lead_brick_8in/2.0 - 6*inch - 3*inch),
                                                lead_behind_radiator_L,
                                                "LeadBehindRadiator",
                                                hall_L,
                                                false,
                                                0);

  // Beam pipe (also redone when radiator geometry was changed)
  if (oldgeo)
  {
    G4Tubs * beampipe_out_S = new G4Tubs("beampipe_out_S", 0, 25*mm, 8*inch, 0, 2.0*M_PI);
    G4LogicalVolume *beampipe_out_L = new G4LogicalVolume(beampipe_out_S, StainlessSteel, "beampipe_out_L");
    G4Tubs * beampipe_in_S = new G4Tubs("beampipe_in_S", 0, 20*mm, 8*inch, 0, 2.0*M_PI);
    G4LogicalVolume *beampipe_in_L = new G4LogicalVolume(beampipe_in_S, Vacuum, "beampipe_in_L");
    G4PVPlacement *lbeampipe_in_P = new G4PVPlacement(0,
                                               G4ThreeVector(0,0,0),
                                               beampipe_in_L,
                                               "BeampipeVacuum",
                                               beampipe_out_L,
                                               false,
                                               0);
    G4PVPlacement *lbeampipe_out_P = new G4PVPlacement(0,
                                               G4ThreeVector(0,0,-8*inch-CuZ/2.0),
                                               beampipe_out_L,
                                               "Beampipe",
                                               hall_L,
                                               false,
                                               0);
  }

  // lead shielding below the radiator
  G4Box *lead_below_radiator_S = new G4Box("lead_below_radiator_S", 5*lead_brick_4in/2.0, 2*lead_brick_2in/2.0, 4*lead_brick_4in/2.0);
  G4LogicalVolume *lead_below_radiator_L = new G4LogicalVolume(lead_below_radiator_S, Pb, "lead_below_radiator_L");
  G4PVPlacement *lead_below_radiator_P = new G4PVPlacement(0,
                                               G4ThreeVector(0,-1.5*lead_brick_4in, - CuZ/2.0 - lead_below_radiator_S->GetZHalfLength() + 2.0*lead_brick_2in),
                                               lead_below_radiator_L,
                                               "LeadBelowRadiator",
                                               hall_L,
                                               false,
                                               0);

  // layer of lead shielding parallel to near detector housing, upbeam
  G4double parallel_layer_det1_thickness = 2.0*inch; // 2" in reality
  G4double parallel_layer_det1_width     = 2.0*lead_brick_8in;

  G4Box *parallel_layer_det1_S = new G4Box("parallel_layer_det1_S", parallel_layer_det1_thickness/2.0, 5*lead_brick_4in/2.0, parallel_layer_det1_width/2.0);
  G4LogicalVolume *parallel_layer_det1_L = new G4LogicalVolume(parallel_layer_det1_S, Pb, "parallel_layer_det1_L");

  //if (useSplitting) biasingOperator->AttachTo(parallel_layer_det1_L);

  G4double parallel_layer_det1_x = det2x - hHPGe_det/2.0 - (lead_filter_width + parallel_layer_det1_thickness/2.0)*sin(relative_angle) + parallel_layer_det1_width/2.0 * cos(relative_angle);
  G4double parallel_layer_det1_y = det1y - det2_h_above_base + 5*lead_brick_4in/2.0;
  G4double parallel_layer_det1_z = det2z - lead_filter_width/2.0 - lead_brick_2in - lead_filter_thickness - (lead_filter_width + parallel_layer_det1_thickness/2.0)*cos(relative_angle) - parallel_layer_det1_width/2.0 * sin(relative_angle);

  G4PVPlacement *parallel_layer_det1_P = new G4PVPlacement(HPGe_rot_1,
                                               G4ThreeVector(parallel_layer_det1_x,parallel_layer_det1_y,parallel_layer_det1_z),
                                               parallel_layer_det1_L,
                                               "ParallelLayerDet1",
                                               hall_L,
                                               false,
                                               0);

  G4PVPlacement *parallel_layer_det1_P2 = new G4PVPlacement(nHPGe_rot_1,
                                               G4ThreeVector(-parallel_layer_det1_x+1.43*inch,parallel_layer_det1_y,parallel_layer_det1_z),
                                               parallel_layer_det1_L,
                                               "ParallelLayerDet12",
                                               hall_L,
                                               false,
                                               0);

  // layer of lead shielding parallel to far detector housing, upbeam
  G4double parallel_layer_det2_thickness = 2.0*lead_brick_2in + lead_filter_thickness;
  G4double parallel_layer_det2_width     = 4.0*lead_brick_4in;
  G4double parallel_layer_det2_height    = 2.0*lead_brick_8in;

  G4Box *parallel_layer_det2_S = new G4Box("parallel_layer_det2_S", parallel_layer_det2_width/2.0, parallel_layer_det2_height/2.0, parallel_layer_det2_thickness/2.0);
  G4LogicalVolume *parallel_layer_det2_L = new G4LogicalVolume(parallel_layer_det2_S, AltPb, "parallel_layer_det2_L");

  G4double parallel_layer_det2_x = det2x - hHPGe_det/2.0 + parallel_layer_det2_width/2.0;
  G4double parallel_layer_det2_y = det2y - det2_h_above_base + parallel_layer_det2_height/2.0;
  G4double parallel_layer_det2_z = det2z - 2.5*inch - parallel_layer_det2_thickness/2.0;

  G4PVPlacement *parallel_layer_det2_P = new G4PVPlacement(0,
                                              G4ThreeVector(parallel_layer_det2_x, parallel_layer_det2_y, parallel_layer_det2_z),
                                              parallel_layer_det2_L,
                                              "ParallelLayerDet2",
                                              hall_L,
                                              false,
                                              0);

  G4PVPlacement *parallel_layer_det2_P2 = new G4PVPlacement(0,
                                              G4ThreeVector(-parallel_layer_det2_x+1.43*inch, parallel_layer_det2_y, parallel_layer_det2_z),
                                              parallel_layer_det2_L,
                                              "ParallelLayerDet22",
                                              hall_L,
                                              false,
                                              0);

  // extra layer of lead shielding by parallel_layer_det2_P
  G4Box *parallel_layer_det2_extra_S = new G4Box("parallel_layer_det2_extra_S",12*inch/2.0, 16*inch/2.0, 2*inch/2.0);
  G4LogicalVolume *parallel_layer_det2_extra_L = new G4LogicalVolume(parallel_layer_det2_extra_S, AltPb, "parallel_layer_det2_extra_L");
  G4PVPlacement *parallel_layer_det2_extra_P = new G4PVPlacement(0,
                                                    G4ThreeVector(parallel_layer_det2_x+2*inch, parallel_layer_det2_y, parallel_layer_det2_z - 3.25*inch),
                                                    parallel_layer_det2_extra_L,
                                                    "ParallelLayerDet2Extra",
                                                    hall_L,
                                                    false,
                                                    0);

  G4PVPlacement *parallel_layer_det2_extra_P2 = new G4PVPlacement(0,
                                                    G4ThreeVector(-parallel_layer_det2_x-2*inch+1.43*inch, parallel_layer_det2_y, parallel_layer_det2_z - 3.25*inch),
                                                    parallel_layer_det2_extra_L,
                                                    "ParallelLayerDet2Extra2",
                                                    hall_L,
                                                    false,
                                                    0);

  // conical lead collimator: a 4x4x8" lead brick with 1/2" and 2" diameter openings
  // can be placed as close as ~3/4" from the radiator face due to the ionization chamber
  G4Box *collimator_lead_S = new G4Box("collimator_lead_S", lead_brick_4in/2.0, lead_brick_4in/2.0, lead_brick_8in/2.0);


  G4double collimator_lead_x = 0.0;
  G4double collimator_lead_y = 0.0;
  collimator_lead_z = CuZ/2+3.0+lead_brick_8in/2.0;//0.75*inch + lead_brick_8in/2.0;

  // conical hole in the lead
  G4double collimator_hole_r1 = 9.86/2.0; // Entry
  G4double collimator_hole_r2 = 26.72/2.0; // Exit

  G4double collimator_hole_h = lead_brick_8in+0.1;
  G4Cons * collimator_hole_S = new G4Cons("collimator_hole_S", 0.0, collimator_hole_r1, 0.0, collimator_hole_r2, collimator_hole_h/2.0, 0, 2.0*M_PI);

  G4SubtractionSolid * collimator_S = new G4SubtractionSolid("collimator_S",collimator_lead_S,collimator_hole_S,0,G4ThreeVector(0,0,0));

  G4LogicalVolume *collimator_lead_L = new G4LogicalVolume(collimator_S, Pb, "collimator_lead_L");

//  G4LogicalVolume *collimator_hole_L = new G4LogicalVolume(collimator_hole_S, Air, "collimator_hole_L");

//  G4PVPlacement *collimator_hole_P = new G4PVPlacement(0,
//                                          G4ThreeVector(),
//                                          collimator_hole_L,
//                                          "CollimatorHole",
//                                          collimator_lead_L,
//                                          false,
//                                          0);

  G4PVPlacement *collimator_lead_P = new G4PVPlacement(0,
                                          G4ThreeVector(collimator_lead_x, collimator_lead_y, collimator_lead_z),
                                          collimator_lead_L,
                                          "CollimatorLead",
                                          hall_L,
                                          false,
                                          0);

  // lead base on which all the shielding sits (just the single bottom layer of bricks)
  G4Box *lead_base_S = new G4Box("lead_base_S", 10*lead_brick_4in/2.0, lead_brick_2in/2.0, 6*lead_brick_8in/2.0);
  G4LogicalVolume *lead_base_L = new G4LogicalVolume(lead_base_S, Pb, "lead_base_L");

  G4double lead_base_x = 2.5*lead_brick_4in + 5*lead_brick_4in;
  G4double lead_base_y = det2y - det2_h_above_base - lead_brick_2in/2.0;
  G4double lead_base_z = lead_behind_radiator_P->GetObjectTranslation().z() + 5*inch + 3*lead_brick_8in;

  G4PVPlacement *lead_base_P = new G4PVPlacement(0,
                                    G4ThreeVector(lead_base_x, lead_base_y, lead_base_z),
                                    lead_base_L,
                                    "LeadBase",
                                    hall_L,
                                    false,
                                    0);

  G4PVPlacement *lead_base_P2 = new G4PVPlacement(0,
                                    G4ThreeVector(-lead_base_x, lead_base_y, lead_base_z),
                                    lead_base_L,
                                    "LeadBase2",
                                    hall_L,
                                    false,
                                    0);

  // tall layer of lead parallel to beamline by near detector
  G4double beamline_lead_tall_thickness = lead_brick_4in;
  G4double beamline_lead_tall_length    = 2.0*lead_brick_8in;
  G4double beamline_lead_tall_height    = 7.0*lead_brick_2in;

  G4Box *beamline_lead_tall_S = new G4Box("beamline_lead_tall_S", beamline_lead_tall_thickness/2.0, beamline_lead_tall_height/2.0, beamline_lead_tall_length/2.0);
  G4LogicalVolume *beamline_lead_tall_L = new G4LogicalVolume(beamline_lead_tall_S, Pb, "beamline_lead_tall_L");

  G4double beamline_lead_tall_x = lead_base_x - 10.0*lead_brick_4in/2.0 + lead_brick_8in;
  G4double beamline_lead_tall_y = lead_base_y + 1.5*lead_brick_2in + beamline_lead_tall_height/2.0;
  G4double beamline_lead_tall_z = lead_base_z - 6.0*lead_brick_8in/2.0 + beamline_lead_tall_length/2.0;

  // the following are currently commented out in order to make room for the giant_lead_shield!

  /*G4PVPlacement *beamline_lead_tall_P = new G4PVPlacement(0,
                                        G4ThreeVector(beamline_lead_tall_x, beamline_lead_tall_y, beamline_lead_tall_z),
                                        beamline_lead_tall_L,
                                        "BeamlineLeadTall",
                                        hall_L,
                                        false,
                                        0);*/

  // short layer of lead parallel to beamline by near detector, excluding the single lead brick standing upright beside it
  G4double beamline_lead_short_thickness = lead_brick_2in;
  G4double beamline_lead_short_length    = lead_brick_4in + lead_brick_8in;
  G4double beamline_lead_short_height    = lead_brick_8in;

  G4Box *beamline_lead_short_S = new G4Box("beamline_lead_short_S", beamline_lead_short_thickness/2.0, beamline_lead_short_height/2.0, beamline_lead_short_length/2.0);
  G4LogicalVolume *beamline_lead_short_L = new G4LogicalVolume(beamline_lead_short_S, Pb, "beamline_lead_short_L");

  G4double beamline_lead_short_x = beamline_lead_tall_x - 1.5*lead_brick_2in;
  G4double beamline_lead_short_y = lead_base_y + 1.5*lead_brick_2in + beamline_lead_short_height/2.0;
  G4double beamline_lead_short_z = lead_base_z - 6.0*lead_brick_8in/2.0 + 2.0*lead_brick_4in + beamline_lead_short_length/2.0;

  /*G4PVPlacement *beamline_lead_short_P = new G4PVPlacement(0,
                                        G4ThreeVector(beamline_lead_short_x, beamline_lead_short_y, beamline_lead_short_z),
                                        beamline_lead_short_L,
                                        "BeamlineLeadShort",
                                        hall_L,
                                        false,
                                        0);*/

  // that one lead brick standing upright by the lead layers parallel to the beamline, since it's pretty close to the radiator
  G4Box *extra_brick_1_S = new G4Box("extra_brick_1_S", lead_brick_4in/2.0, lead_brick_8in/2.0, lead_brick_2in/2.0);
  G4LogicalVolume *extra_brick_1_L = new G4LogicalVolume(extra_brick_1_S, Pb, "extra_brick_1_L");

  /*G4PVPlacement *extra_brick_1_P = new G4PVPlacement(0,
                                        G4ThreeVector(beamline_lead_short_x - 1*inch, beamline_lead_short_y, beamline_lead_short_z - 7*inch),
                                        extra_brick_1_L,
                                        "ExtraBrick1",
                                        hall_L,
                                        false,
                                        0);*/

  // extra base layer of bricks under beamline_lead_{short,tall}
  G4double lead_base_extra_length = 2*lead_brick_8in + lead_brick_4in;
  G4double lead_base_extra_width  = 3.4*lead_brick_4in; // hard to judge from photos, but whatever works
  G4double lead_base_extra_height = lead_brick_2in;
  G4Box *lead_base_extra_S = new G4Box("lead_base_extra_S", lead_base_extra_length/2.0, lead_base_extra_height/2.0, lead_base_extra_width/2.0);
  G4LogicalVolume *lead_base_extra_L = new G4LogicalVolume(lead_base_extra_S, Pb, "lead_base_extra_L");

  G4double lead_base_extra_x = lead_base_x - 10.0*lead_brick_4in/2.0 + lead_base_extra_length/2.0;
  G4double lead_base_extra_y = lead_base_y + lead_brick_2in;
  G4double lead_base_extra_z = lead_base_z - 6.0*lead_brick_8in/2.0 + lead_base_extra_width/2.0;

  G4PVPlacement *lead_base_extra_P = new G4PVPlacement(0,
                                          G4ThreeVector(lead_base_extra_x, lead_base_extra_y, lead_base_extra_z),
                                          lead_base_extra_L,
                                          "LeadBaseExtra",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_base_extra_P2 = new G4PVPlacement(0,
                                          G4ThreeVector(-lead_base_extra_x, lead_base_extra_y, lead_base_extra_z),
                                          lead_base_extra_L,
                                          "LeadBaseExtra2",
                                          hall_L,
                                          false,
                                          0);

  // roughly wedge-shaped base that must be buried under all the rest of the lead (since I can't see it directly)
  // x,y,z don't really match up to x,y,z directions due to rotation
  G4double lead_wedge_dx1 = 0.01*mm; // short x dimension; exactly 0 is prohibited
  G4double lead_wedge_dx2 = 15.0*inch; // long x dimension (best coverage)
  G4double lead_wedge_dz  = lead_brick_2in;
  G4double lead_wedge_dy  = lead_wedge_dx2 * tan(relative_angle);

  G4Trap *lead_wedge_S = new G4Trap("lead_wedge_S", lead_wedge_dz, lead_wedge_dy, lead_wedge_dx2, lead_wedge_dx1);
  G4LogicalVolume *lead_wedge_L = new G4LogicalVolume(lead_wedge_S, Pb, "lead_wedge_L");

  G4RotationMatrix *wedge_rot = new G4RotationMatrix();
  wedge_rot->rotateX(3*pi/2.0);

  G4RotationMatrix *nwedge_rot = new G4RotationMatrix();
  nwedge_rot->rotateX(pi/2.0);
  nwedge_rot->rotateZ(pi);

  G4double lead_wedge_x = lead_base_extra_x - lead_base_extra_length/2.0 + lead_wedge_dx2/2.0 - 9.5*cm; // fudge factor for weird x coord
  G4double lead_wedge_y = lead_base_y + lead_brick_2in;
  G4double lead_wedge_z = lead_base_extra_z + lead_base_extra_width/2.0 + lead_wedge_dy/2.0;

  G4PVPlacement *lead_wedge_P = new G4PVPlacement(wedge_rot,
                                      G4ThreeVector(lead_wedge_x, lead_wedge_y, lead_wedge_z),
                                      lead_wedge_L,
                                      "LeadWedge",
                                      hall_L,
                                      false,
                                      0);

  G4PVPlacement *lead_wedge_P2 = new G4PVPlacement(nwedge_rot,
                                      G4ThreeVector(-lead_wedge_x+1.43*inch, lead_wedge_y, lead_wedge_z),
                                      lead_wedge_L,
                                      "LeadWedge2",
                                      hall_L,
                                      false,
                                      0);

  // (Presumably) stainless steel housings the detectors sit in. Ignore the rounded corners, despite whatever Apple says.
  // Because I've already placed the HPGe detectors and they're a vital part of my coordinate system
  // (and the housings are open-ended on either side), hack together something using subtraction volumes
  G4double ss_housing_thickness = 0.5*cm;
  G4double ss_housing_length = 2.0*lead_brick_8in;
  G4double ss_housing_width  = 2.5*lead_brick_2in;
  G4double ss_housing_height = ss_housing_width;

  G4Box *ss_outer_S = new G4Box("ss_outer_S", ss_housing_width/2.0, ss_housing_height/2.0, ss_housing_length/2.0);

  G4Box *ss_inner_S = new G4Box("ss_inner_S",
    (ss_housing_width  - ss_housing_thickness)/2.0,
    (ss_housing_height - ss_housing_thickness)/2.0,
    (ss_housing_length)/2.0+40);

  G4SubtractionSolid *ss_housing_S = new G4SubtractionSolid("ss_housing_S", ss_outer_S, ss_inner_S);
  G4LogicalVolume *ss_housing_L1 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L1");
  G4LogicalVolume *ss_housing_L2 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L2");
  G4LogicalVolume *ss_housing_L3 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L3");
  G4LogicalVolume *ss_housing_L4 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L4");
  G4LogicalVolume *ss_housing_L5 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L5");
  G4LogicalVolume *ss_housing_L6 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L6");
  G4LogicalVolume *ss_housing_L7 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L7");
  G4LogicalVolume *ss_housing_L8 = new G4LogicalVolume(ss_housing_S, StainlessSteel, "ss_housing_L8");

  G4double ss_housing_det1_x = det1x - hHPGe_det/2.0 * cos(relative_angle) + ss_housing_length/2.0 * cos(relative_angle);
  G4double ss_housing_det1_y = det1y;
  G4double ss_housing_det1_z = det1z + hHPGe_det/2.0 * sin(relative_angle) - ss_housing_length/2.0 * sin(relative_angle);

  if (incfilters)
  {
    G4PVPlacement *ss_housing_1_P = new G4PVPlacement(HPGe_rot_1,
                                          G4ThreeVector(ss_housing_det1_x, ss_housing_det1_y, ss_housing_det1_z),
                                          ss_housing_L1,
                                          "StainlessSteelHousing1",
                                          hall_L,
                                          false,
                                          0);

    G4PVPlacement *ss_housing_1_P2 = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-ss_housing_det1_x+1.43*inch, ss_housing_det1_y, ss_housing_det1_z),
                                          ss_housing_L2,
                                          "StainlessSteelHousing12",
                                          hall_L,
                                          false,
                                          0);
  }

  G4double ss_housing_det2_x = det2x - hHPGe_det/2.0 + ss_housing_length/2.0;
  G4double ss_housing_det2_y = det2y;
  G4double ss_housing_det2_z = det2z;

  if (incfilters)
  {
    G4PVPlacement *ss_housing_2_P = new G4PVPlacement(HPGe_rot_2,
                                          G4ThreeVector(ss_housing_det2_x, ss_housing_det2_y, ss_housing_det2_z),
                                          ss_housing_L3,
                                          "StainlessSteelHousing2",
                                          hall_L,
                                          false,
                                          0);

    G4PVPlacement *ss_housing_2_P2 = new G4PVPlacement(nHPGe_rot_2,
                                          G4ThreeVector(-ss_housing_det2_x+1.43*inch, ss_housing_det2_y, ss_housing_det2_z),
                                          ss_housing_L4,
                                          "StainlessSteelHousing22",
                                          hall_L,
                                          false,
                                          0);

  }


  G4double ss_housing_det3_x = ss_housing_det1_x;
  G4double ss_housing_det3_y = ss_housing_det1_y + ss_housing_height;
  G4double ss_housing_det3_z = ss_housing_det1_z;

  if (incfilters)
  {
    G4PVPlacement *ss_housing_3_P = new G4PVPlacement(HPGe_rot_1,
                                          G4ThreeVector(ss_housing_det3_x, ss_housing_det3_y, ss_housing_det3_z),
                                          ss_housing_L5,
                                          "StainlessSteelHousing3",
                                          hall_L,
                                          false,
                                          0);

    G4PVPlacement *ss_housing_3_P2 = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-ss_housing_det3_x+1.43*inch, ss_housing_det3_y, ss_housing_det3_z),
                                          ss_housing_L6,
                                          "StainlessSteelHousing32",
                                          hall_L,
                                          false,
                                          0);

    // Missing detector phantom brick
    G4Box * detphantombrick_S = new G4Box("detphantombrick_S", lead_brick_4in/2.0, lead_brick_2in/2.0, lead_brick_8in/2.0);
    G4LogicalVolume * detphantombrick_L = new G4LogicalVolume(detphantombrick_S, Pb, "detphantombrick_L");
    G4PVPlacement * detphantombrick_P = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-ss_housing_det3_x + 4.0*25.4*cos(35.0/360.0*2.0*M_PI) +1.43*inch, ss_housing_det3_y-35.6, ss_housing_det3_z + + 4.0*25.4*sin(35.0/360.0*2.0*M_PI)),
                                          detphantombrick_L,
                                          "detphantombrick",
                                          hall_L,
                                          false,
                                          0);

  }

  G4double ss_housing_det4_x = ss_housing_det2_x;
  G4double ss_housing_det4_y = ss_housing_det2_y + ss_housing_height;
  G4double ss_housing_det4_z = ss_housing_det2_z;

  if (incfilters)
  {
    G4PVPlacement *ss_housing_4_P = new G4PVPlacement(HPGe_rot_2,
                                          G4ThreeVector(ss_housing_det4_x, ss_housing_det4_y, ss_housing_det4_z),
                                          ss_housing_L7,
                                          "StainlessSteelHousing4",
                                          hall_L,
                                          false,
                                          0);

    G4PVPlacement *ss_housing_4_P2 = new G4PVPlacement(nHPGe_rot_2,
                                          G4ThreeVector(-ss_housing_det4_x+1.43*inch, ss_housing_det4_y, ss_housing_det4_z),
                                          ss_housing_L8,
                                          "StainlessSteelHousing42",
                                          hall_L,
                                          false,
                                          0);
  }

  G4VisAttributes *HousingVisAtt = new G4VisAttributes(G4Color(50./256, 70./256, 245./256));
  HousingVisAtt->SetVisibility(true);
  ss_housing_L1->SetVisAttributes(HousingVisAtt);
  ss_housing_L2->SetVisAttributes(HousingVisAtt);

  // 12x16" lead shield downbeam of detector 2
  G4Box *downbeam_lead_det2_S = new G4Box("downbeam_lead_det2_S", 2.0*lead_brick_8in/2.0, 3.0*lead_brick_4in/2.0, 2.0*lead_brick_2in/2.0);
  G4LogicalVolume *downbeam_lead_det2_L = new G4LogicalVolume(downbeam_lead_det2_S, Pb, "downbeam_lead_det2_L");

  G4double downbeam_lead_det2_x = ss_housing_det2_x;
  G4double downbeam_lead_det2_y = lead_base_y + lead_brick_2in/2.0 + 1.5*lead_brick_4in;
  G4double downbeam_lead_det2_z = ss_housing_det2_z + ss_housing_width/2.0 + lead_brick_2in;

  G4PVPlacement *downbeam_lead_det2_P = new G4PVPlacement(0,
                                              G4ThreeVector(downbeam_lead_det2_x, downbeam_lead_det2_y, downbeam_lead_det2_z),
                                              downbeam_lead_det2_L,
                                              "DownbeamLead2",
                                              hall_L,
                                              false,
                                              0);

  G4PVPlacement *downbeam_lead_det2_P2 = new G4PVPlacement(0,
                                              G4ThreeVector(-downbeam_lead_det2_x+1.43*inch, downbeam_lead_det2_y, downbeam_lead_det2_z),
                                              downbeam_lead_det2_L,
                                              "DownbeamLead22",
                                              hall_L,
                                              false,
                                              0);

  // 1" thick lead shims that sit above and below the blue stainless housings
  G4Box *lead_shim_S = new G4Box("lead_shim_S", ss_housing_width/2.0, 1.0*inch/2.0, ss_housing_length/2.0);
  G4LogicalVolume *lead_shim_L = new G4LogicalVolume(lead_shim_S, Pb, "lead_shim_L");
  G4LogicalVolume *lead_shim_L2 = new G4LogicalVolume(lead_shim_S, Pb, "lead_shim_L2");

  G4PVPlacement *lead_shim_det2_P = new G4PVPlacement(HPGe_rot_2,
                                          G4ThreeVector(ss_housing_det2_x, ss_housing_det2_y - 5.0*inch/2.0 - 1.0*inch/2.0, ss_housing_det2_z),
                                          lead_shim_L2,
                                          "LeadShim2",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_shim_det2_P2 = new G4PVPlacement(nHPGe_rot_2,
                                          G4ThreeVector(-ss_housing_det2_x+1.43*inch, ss_housing_det2_y - 5.0*inch/2.0 - 1.0*inch/2.0, ss_housing_det2_z),
                                          lead_shim_L2,
                                          "LeadShim22",
                                          hall_L,
                                          false,
                                          0);

  G4double lead_shim_det4_y = ss_housing_det4_y + 5.0*inch/2.0 + 1.0*inch/2.0;
  G4PVPlacement *lead_shim_det4_P = new G4PVPlacement(HPGe_rot_2,
                                          G4ThreeVector(ss_housing_det4_x, lead_shim_det4_y, ss_housing_det4_z),
                                          lead_shim_L2,
                                          "LeadShim4",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_shim_det4_P2 = new G4PVPlacement(nHPGe_rot_2,
                                          G4ThreeVector(-ss_housing_det4_x+1.43*inch, lead_shim_det4_y, ss_housing_det4_z),
                                          lead_shim_L2,
                                          "LeadShim42",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_shim_det1_P = new G4PVPlacement(HPGe_rot_1,
                                          G4ThreeVector(ss_housing_det1_x, ss_housing_det1_y - 5.0*inch/2.0 - 1.0*inch/2.0, ss_housing_det1_z),
                                          lead_shim_L,
                                          "LeadShim1",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_shim_det1_P2 = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-ss_housing_det1_x+1.43*inch, ss_housing_det1_y - 5.0*inch/2.0 - 1.0*inch/2.0, ss_housing_det1_z),
                                          lead_shim_L,
                                          "LeadShim12",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_shim_det3_P = new G4PVPlacement(HPGe_rot_1,
                                          G4ThreeVector(ss_housing_det3_x, ss_housing_det3_y + 5.0*inch/2.0 + 1.0*inch/2.0, ss_housing_det3_z),
                                          lead_shim_L,
                                          "LeadShim3",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_shim_det3_P2 = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-ss_housing_det3_x+1.43*inch, ss_housing_det3_y + 5.0*inch/2.0 + 1.0*inch/2.0, ss_housing_det3_z),
                                          lead_shim_L,
                                          "LeadShim32",
                                          hall_L,
                                          false,
                                          0);

  // lead roof on detector 2 housing (downbeam)
  G4Box *lead_roof_det2_S = new G4Box("lead_roof_det2_S", 2.0*lead_brick_8in/2.0, lead_brick_4in/2.0, 4.0*lead_brick_2in/2.0);
  G4LogicalVolume *lead_roof_det2_L = new G4LogicalVolume(lead_roof_det2_S, Pb, "lead_roof_det2_L");

  G4PVPlacement *lead_roof_det2_P = new G4PVPlacement(0,
                                          G4ThreeVector(ss_housing_det2_x, lead_shim_det4_y + 0.5*inch + lead_brick_4in/2.0, ss_housing_det2_z - 2.5*inch + 2.0*lead_brick_2in),
                                          lead_roof_det2_L,
                                          "LeadRoof2",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_roof_det2_P2 = new G4PVPlacement(0,
                                          G4ThreeVector(-ss_housing_det2_x+1.43*inch, lead_shim_det4_y + 0.5*inch + lead_brick_4in/2.0, ss_housing_det2_z - 2.5*inch + 2.0*lead_brick_2in),
                                          lead_roof_det2_L,
                                          "LeadRoof22",
                                          hall_L,
                                          false,
                                          0);


  // detector 1 shielding, upbeam; sits between parallel layer det1 and the stainless steel housing
  G4double upbeam_lead_det1_length = 2.0*lead_brick_8in;
  G4double upbeam_lead_det1_height = lead_brick_8in + lead_brick_2in + lead_brick_2in;
  G4double upbeam_lead_det1_thickness = lead_brick_2in - 0.5*inch;

  G4Box *upbeam_lead_det1_S = new G4Box("upbeam_lead_det1_S", upbeam_lead_det1_thickness/2.0, upbeam_lead_det1_height/2.0, upbeam_lead_det1_length/2.0);
  G4LogicalVolume *upbeam_lead_det1_L = new G4LogicalVolume(upbeam_lead_det1_S, Pb, "upbeam_lead_det1_L");

  //if (useSplitting) biasingOperator->AttachTo(upbeam_lead_det1_L);

  G4double upbeam_lead_det1_x = parallel_layer_det1_x + 0.5*(parallel_layer_det1_thickness + upbeam_lead_det1_thickness) * sin(relative_angle);
  G4double upbeam_lead_det1_y = lead_base_y + lead_brick_2in/2.0 + upbeam_lead_det1_height/2.0;
  G4double upbeam_lead_det1_z = parallel_layer_det1_z + 0.5*(parallel_layer_det1_thickness + upbeam_lead_det1_thickness) * cos(relative_angle);

  G4VPhysicalVolume *upbeam_lead_det1_P = new G4PVPlacement(HPGe_rot_1,
                                                G4ThreeVector(upbeam_lead_det1_x, upbeam_lead_det1_y, upbeam_lead_det1_z),
                                                upbeam_lead_det1_L,
                                                "UpbeamLeadDet1",
                                                hall_L,
                                                false,
                                                0);

  G4VPhysicalVolume *upbeam_lead_det1_P2 = new G4PVPlacement(nHPGe_rot_1,
                                                G4ThreeVector(-upbeam_lead_det1_x+1.43*inch, upbeam_lead_det1_y, upbeam_lead_det1_z),
                                                upbeam_lead_det1_L,
                                                "UpbeamLeadDet12",
                                                hall_L,
                                                false,
                                                0);

  // detector 1 shielding, downbeam; sits between the upbeam-most lead layer of det2 and the downbeam side of det2 stainless housing
  G4double corner_cutoff = 1.0*inch; // cut off a bit to prevent ovelaps
  G4double downbeam_lead_det1_thickness = upbeam_lead_det1_thickness;
  G4double downbeam_lead_det1_length = upbeam_lead_det1_length - corner_cutoff;
  G4double downbeam_lead_det1_height = upbeam_lead_det1_height;

  G4Box *downbeam_lead_det1_S = new G4Box("upbeam_lead_det1_S", downbeam_lead_det1_thickness/2.0, downbeam_lead_det1_height/2.0, downbeam_lead_det1_length/2.0);

  G4LogicalVolume *downbeam_lead_det1_L = new G4LogicalVolume(downbeam_lead_det1_S, AltPb, "downbeam_lead_det1_L");

  G4double downbeam_lead_det1_x = ss_housing_det1_x + 0.5*(ss_housing_width + downbeam_lead_det1_thickness) * sin(relative_angle) + corner_cutoff/2.0 * cos(relative_angle);
  G4double downbeam_lead_det1_y = upbeam_lead_det1_y;
  G4double downbeam_lead_det1_z = ss_housing_det1_z + 0.5*(ss_housing_width + downbeam_lead_det1_thickness) * cos(relative_angle) - corner_cutoff/2.0 * sin(relative_angle);

  G4VPhysicalVolume *downbeam_lead_det1_P = new G4PVPlacement(HPGe_rot_1,
                                                  G4ThreeVector(downbeam_lead_det1_x, downbeam_lead_det1_y, downbeam_lead_det1_z),
                                                  downbeam_lead_det1_L,
                                                  "DownbeamLeadDet1",
                                                  hall_L,
                                                  false,
                                                  0);

  G4VPhysicalVolume *downbeam_lead_det1_P2 = new G4PVPlacement(nHPGe_rot_1,
                                                  G4ThreeVector(-downbeam_lead_det1_x+1.43*inch, downbeam_lead_det1_y, downbeam_lead_det1_z),
                                                  downbeam_lead_det1_L,
                                                  "DownbeamLeadDet12",
                                                  hall_L,
                                                  false,
                                                  0);


  //G4cout<<"\n\n"<<G4ThreeVector(downbeam_lead_det1_x, downbeam_lead_det1_y, downbeam_lead_det1_z)<<" "<<G4ThreeVector(upbeam_lead_det1_x, upbeam_lead_det1_y, upbeam_lead_det1_z)<<"\n\n";


  // lead roof on detector 1 housing; there is at least 8x4x16" of lead everywhere, and more in some places
  // subtract off 1 inch from width to prevent overlaps
  G4Box *lead_roof_det1_S = new G4Box("lead_roof_det1_S", (lead_brick_8in-1*inch)/2.0, lead_brick_4in/2.0, 2.0*lead_brick_8in/2.0);
  G4LogicalVolume *lead_roof_det1_L = new G4LogicalVolume(lead_roof_det1_S, Pb, "lead_roof_det1_L");

  G4double lead_roof_det1_x = ss_housing_det1_x;
  G4double lead_roof_det1_y = downbeam_lead_det1_y + downbeam_lead_det1_height/2.0 + lead_brick_4in/2.0;
  G4double lead_roof_det1_z = ss_housing_det1_z;

  G4PVPlacement *lead_roof_det1_P = new G4PVPlacement(HPGe_rot_1,
                                          G4ThreeVector(lead_roof_det1_x, lead_roof_det1_y, lead_roof_det1_z),
                                          lead_roof_det1_L,
                                          "LeadRoof1",
                                          hall_L,
                                          false,
                                          0);

  G4PVPlacement *lead_roof_det1_P2 = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-lead_roof_det1_x+1.43*inch, lead_roof_det1_y, lead_roof_det1_z),
                                          lead_roof_det1_L,
                                          "LeadRoof12",
                                          hall_L,
                                          false,
                                          0);

  // ------------------------------------------------------------------------------------
  // Optional lead shielding
  //

  // giant lead shield in front of upbeam detector
  // build mother volume out of vacuum, then fill with lead layers for splitting algorithm
  // also geant is silly, so pull some tricks with the kXAxis and rotation matrix to make it work
  G4int nLayers = 4;
  G4double layerThickness = 2.0*inch;
  G4double shieldHeight = 2.0*lead_brick_8in;
  G4double shieldLength = 2.0*lead_brick_8in;
  G4double totalThickness = nLayers*layerThickness;

  G4double shieldX = parallel_layer_det1_x - sin(relative_angle) * (parallel_layer_det1_thickness + totalThickness)/2.0;
  G4double shieldY = lead_wedge_y + lead_brick_2in/2.0 + shieldHeight/2.0;
  G4double shieldZ = parallel_layer_det1_z - cos(relative_angle) * (parallel_layer_det1_thickness + totalThickness)/2.0;

  G4Box *giant_lead_shield_S = new G4Box("giant_lead_shield_S", totalThickness/2.0, shieldHeight/2.0, shieldLength/2.0);
  G4LogicalVolume *giant_lead_shield_L = new G4LogicalVolume(giant_lead_shield_S, Vacuum, "giant_lead_shield_L");
  G4PVPlacement *giant_lead_shield_P = new G4PVPlacement(HPGe_rot_1,
                                            G4ThreeVector(shieldX, shieldY, shieldZ),
                                            giant_lead_shield_L,
                                            "GiantLeadShield",
                                            hall_L,
                                            false,
                                            0);

  G4PVPlacement *giant_lead_shield_P2 = new G4PVPlacement(nHPGe_rot_1,
                                            G4ThreeVector(-shieldX+1.43*inch, shieldY, shieldZ),
                                            giant_lead_shield_L,
                                            "GiantLeadShield2",
                                            hall_L,
                                            false,
                                            0);


  G4Box *lead_layer_S = new G4Box("lead_layer_S", layerThickness/2.0, shieldHeight/2.0, shieldLength/2.0);
  G4LogicalVolume *lead_layer_L = new G4LogicalVolume(lead_layer_S, Pb, "lead_layer_L");

  //if (useSplitting) biasingOperator->AttachTo(lead_layer_L);

  G4PVReplica *lead_layer_P = new G4PVReplica("LeadLayer",
                                    lead_layer_L,
                                    giant_lead_shield_L,
                                    kXAxis,
                                    nLayers,
                                    layerThickness,
                                    0);

  //G4VisAttributes *GiantLeadShieldVisAtt = new G4VisAttributes(G4Color(255/256.,50/256.,50/256.));
  //giant_lead_shield_L->SetVisAttributes(GiantLeadShieldVisAtt);
  //lead_layer_L->SetVisAttributes(GiantLeadShieldVisAtt);


  // it might be useful to put some lead in the beamline (post-target) to clean up some scatters off the far wall
  G4Box *lead_beamdump_S = new G4Box("lead_beamdump_S", 2.0*lead_brick_8in/2.0, 8.0*lead_brick_2in/2.0, lead_brick_4in/2.0);
  G4LogicalVolume *lead_beamdump_L = new G4LogicalVolume(lead_beamdump_S, Pb, "lead_beamdump_L");

  /*G4PVPlacement *lead_beamdump_P = new G4PVPlacement(0,
                                        G4ThreeVector(0,0,2.0*meter),
                                        lead_beamdump_L,
                                        "LeadBeamDump",
                                        hall_L,
                                        false,
                                        0);*/


  // additional possible shielding around the radiator, because the collimator doesn't catch everything
  // four lead bricks (say) stacked vertically (so 8" height) with their long sides parallel to the beam line
  // for now, put them flush with the shielding behind the radiator
  G4Box *collimator_side_shield_S = new G4Box("collimator_side_shield_S", lead_brick_4in/2.0, 4*lead_brick_2in/2.0, lead_brick_8in/2.0);
  G4Box *collimator_side_shield1_S = new G4Box("collimator_side_shield_S", lead_brick_4in/2.0, lead_brick_8in/2.0, lead_brick_2in/2.0);

  G4LogicalVolume *collimator_shield_right_L = new G4LogicalVolume(collimator_side_shield_S, Pb, "collimator_shield_right_L");
  G4LogicalVolume *collimator_shield_left_L = new G4LogicalVolume(collimator_side_shield_S, Pb, "collimator_shield_left_L");
  G4LogicalVolume *collimator_shield1_right_L = new G4LogicalVolume(collimator_side_shield1_S, Pb, "collimator_shield_right1_L");
  G4LogicalVolume *collimator_shield1_left_L = new G4LogicalVolume(collimator_side_shield1_S, Pb, "collimator_shield_left1_L");

  //if (useSplitting) biasingOperator->AttachTo(collimator_shield_right_L);

  G4double collimator_side_shield_xr =  lead_brick_4in/2.0 + lead_brick_4in/2.0;
  G4double collimator_side_shield_xl = -lead_brick_4in/2.0 - lead_brick_4in/2.0;
  G4double collimator_side_shield_y  = 0.0;
  G4double collimator_side_shield_z = collimator_lead_z; // lead_behind_radiator_P->GetObjectTranslation().z() + 1.5*lead_brick_4in/2.0 + lead_brick_8in/2.0;

  G4PVPlacement *collimator_shield_right_P = new G4PVPlacement(0,
                                                  G4ThreeVector(collimator_side_shield_xr, collimator_side_shield_y, collimator_side_shield_z),
                                                  collimator_shield_right_L,
                                                  "CollimatorShieldRight",
                                                  hall_L,
                                                  false,
                                                  0);

  G4PVPlacement *collimator_shield_left_P = new G4PVPlacement(0,
                                                  G4ThreeVector(collimator_side_shield_xl, collimator_side_shield_y, collimator_side_shield_z),
                                                  collimator_shield_left_L,
                                                  "CollimatorShieldLeft",
                                                  hall_L,
                                                  false,
                                                  0);

  G4PVPlacement *collimator_shield1_right_P = new G4PVPlacement(0,
                                                  G4ThreeVector(collimator_side_shield_xr+15*mm, collimator_side_shield_y, collimator_side_shield_z-5.0*inch),
                                                  collimator_shield1_right_L,
                                                  "CollimatorShieldRight1",
                                                  hall_L,
                                                  false,
                                                  0);

  G4PVPlacement *collimator_shield1_left_P = new G4PVPlacement(0,
                                                  G4ThreeVector(collimator_side_shield_xl-15*mm, collimator_side_shield_y, collimator_side_shield_z-5.0*inch),
                                                  collimator_shield1_left_L,
                                                  "CollimatorShieldLeft1",
                                                  hall_L,
                                                  false,
                                                  0);

  // lead layer above the collimator, sitting on top of the left/right shields: two stacks of two lead bricks
  G4Box *collimator_shield_top_S = new G4Box("collimator_shield_top_S", lead_brick_8in/2.0, 2*lead_brick_2in/2.0, 2*lead_brick_4in/2.0);
  G4LogicalVolume *collimator_shield_top_L = new G4LogicalVolume(collimator_shield_top_S, Pb, "collimator_shield_top_L");

  G4double collimator_shield_top_x = 0;
  G4double collimator_shield_top_y = lead_brick_4in + lead_brick_2in;
  G4double collimator_shield_top_z = lead_behind_radiator_P->GetObjectTranslation().z() + 1.5*lead_brick_4in/2.0 + lead_brick_8in/2.0;

  G4PVPlacement *collimator_shield_top_P = new G4PVPlacement(0,
                                                G4ThreeVector(collimator_shield_top_x, collimator_shield_top_y, collimator_shield_top_z),
                                                collimator_shield_top_L,
                                                "CollimatorShieldTop",
                                                hall_L,
                                                false,
                                                0);

  // slide another lead brick under the collimator, because we'll need one to support the collimator anyway
  G4Box *collimator_shield_bot_S = new G4Box("collimator_shield_bot_S", lead_brick_4in/2.0, lead_brick_2in/2.0*3.0/4.0, lead_brick_8in/2.0);
  G4LogicalVolume *collimator_shield_bot_L = new G4LogicalVolume(collimator_shield_bot_S, Pb, "collimator_shield_bot_L");

  G4PVPlacement *collimator_shield_bot_P = new G4PVPlacement(0,
                                                G4ThreeVector(0, -3.0*inch, collimator_lead_z),
                                                collimator_shield_bot_L,
                                                "CollimatorShieldBottom",
                                                hall_L,
                                                false,
                                                0);

  G4VisAttributes *testVisAtt = new G4VisAttributes();
  testVisAtt->SetForceSolid(true);
  //collimator_shield_bot_L->SetVisAttributes(testVisAtt);


  if (incfilters)
  {
    G4PVPlacement *extra_filter_P = new G4PVPlacement(HPGe_rot_1,
                                          G4ThreeVector(extra_filter_x, extra_filter_y, extra_filter_z),
                                          extra_filter_L,
                                          "ExtraFilter1",
                                          hall_L,
                                          false,
                                          0);

    G4PVPlacement *extra_filter_P2 = new G4PVPlacement(nHPGe_rot_1,
                                          G4ThreeVector(-extra_filter_x+1.43*inch, extra_filter_y, extra_filter_z),
                                          extra_filter_L,
                                          "ExtraFilter12",
                                          hall_L,
                                          false,
                                          0);
  }


  //G4cout<<"\n\nExtra filter pos: "<<G4ThreeVector(extra_filter_x, extra_filter_y, extra_filter_z)<<"\n\n";


  // lead cylinder behind detector 1
  G4double lead_cylinder_r = 2*inch;
  G4double lead_cylinder_h = 4*inch;

  G4Tubs *lead_cylinder_S = new G4Tubs("lead_cylinder_S", 0, lead_cylinder_r, lead_cylinder_h/2.0, 0, 2.0*M_PI);
  G4LogicalVolume *lead_cylinder_L = new G4LogicalVolume(lead_cylinder_S, Pb, "lead_cylinder_L");

  G4double lead_cylinder_x = det1x_actual + (hHPGe_det + lead_cylinder_h + 180)/2.0 * cos(relative_angle);
  G4double lead_cylinder_y = det1y_actual;
  G4double lead_cylinder_z = det1z_actual - (hHPGe_det + lead_cylinder_h + 180)/2.0 * sin(relative_angle);

//  G4PVPlacement *lead_cylinder_P = new G4PVPlacement(HPGe_rot_1,
//                                        G4ThreeVector(lead_cylinder_x, lead_cylinder_y, lead_cylinder_z),
//                                        lead_cylinder_L,
//                                        "LeadCylinder1",
//                                        hall_L,
//                                        false,
//                                        0);

  // extra shielding below det 1
  G4double lead_bottom_height = lead_brick_4in;
  G4double lead_bottom_length = lead_brick_8in;
  G4double lead_bottom_width  = lead_brick_8in;
  G4Box *lead_bottom_S = new G4Box("lead_bottom_S", lead_bottom_width/2.0, lead_bottom_height/2.0, lead_bottom_length/2.0);
  G4LogicalVolume *lead_bottom_L = new G4LogicalVolume(lead_bottom_S, Pb, "lead_bottom_L");

  G4double lead_bottom_x = det1x_actual;
  G4double lead_bottom_y = lead_base_y - lead_brick_2in/2.0 - lead_bottom_height/2.0;
  G4double lead_bottom_z = det1z_actual;

  /*G4PVPlacement *lead_bottom_P = new G4PVPlacement(0,
                                      G4ThreeVector(lead_bottom_x, lead_bottom_y, lead_bottom_z),
                                      lead_bottom_L,
                                      "LeadBottom1",
                                      hall_L,
                                      false,
                                      0);*/

  // corner post
  G4double corner_post_height = 2.0*lead_brick_8in;
  G4double corner_post_length = lead_brick_4in;
  G4double corner_post_thickness = lead_brick_2in;

  G4Box *corner_post_S = new G4Box("corner_post_S", corner_post_thickness/2.0, corner_post_height/2.0, corner_post_length/2.0);
  G4LogicalVolume *corner_post_L = new G4LogicalVolume(corner_post_S, Pb, "corner_post_L");

  G4double corner_post_x = parallel_layer_det1_x - (parallel_layer_det1_width + corner_post_length)/2.0 * cos(relative_angle);
  G4double corner_post_y = lead_base_y + lead_brick_2in/2.0 + corner_post_height/2.0;
  G4double corner_post_z = parallel_layer_det1_z + (parallel_layer_det1_width + corner_post_length)/2.0 * sin(relative_angle);

  G4PVPlacement *corner_post_P = new G4PVPlacement(HPGe_rot_1,
                                        G4ThreeVector(corner_post_x, corner_post_y, corner_post_z),
                                        corner_post_L,
                                        "CornerShield1",
                                        hall_L,
                                        false,
                                        0);

  G4PVPlacement *corner_post_P2 = new G4PVPlacement(nHPGe_rot_1,
                                        G4ThreeVector(-corner_post_x+1.43*inch, corner_post_y, corner_post_z),
                                        corner_post_L,
                                        "CornerShield2",
                                        hall_L,
                                        false,
                                        0);


  corner_post_L->SetVisAttributes(testVisAtt);


  // ------------------------------------------------------------------------------------
  // Uranium (or Al) plate in the beam
  //
  G4Material * platemat = U;
  double plateext = 1.0;
  double platethick = 4.0/32.0*inch;
  std::cout<<"\n\n"<<plateconf<<"\n\n";
  if (plateconf==0 || plateconf==3)
  {
    platemat = U;
    platethick = 4.0/32.0*inch;
    plateext = 1.0;
  }
  else if (plateconf==1)
  {
    platemat = U;
    platethick = 1.0/32.0*inch;
    plateext = 1.0;
  }
  else if (plateconf==2)
  {
    platemat = Al;
    platethick = 2.5*inch;
    plateext = 2.0;
  }

  G4double U_plate_length = 6.0*inch*plateext;
  G4double U_plate_width  = 6.0*inch*plateext; // approximately
  G4double U_plate_thickness = platethick; // thinnest DU plate we have is 1/32"
  //G4cout<<"Platethick: "<<platethick<<"\n\n";

  G4Box *U_plate_S = new G4Box("U_plate_S", U_plate_width/2.0, U_plate_length/2.0, U_plate_thickness/2.0);
  G4LogicalVolume *U_plate_L = new G4LogicalVolume(U_plate_S, platemat, "U_plate_L");

  U_plate_z = 998.9; // det1z_actual + det1x_actual * tan(relative_angle);

  if (plateconf<4)
  {
    std::cout<<"\n\nPlacing U plate at "<<U_plate_z<<"\n\n";
    // exit(0);
    G4PVPlacement *U_plate_P = new G4PVPlacement(0,
                                    G4ThreeVector(0, 0, U_plate_z),
                                    U_plate_L,
                                    "UraniumPlate",
                                    hall_L,
                                    false,
                                    0);
  }

  if (plateconf==3)
  {
    platemat = Al;
    platethick = 2.5*inch;
    plateext = 2.0;

    U_plate_length = 6.0*inch*plateext;
    U_plate_width  = 6.0*inch*plateext; // approximately
    U_plate_thickness = platethick; // thinnest DU plate we have is 1/32"

    G4Box *Al_plate_S = new G4Box("Al_plate_S", U_plate_width/2.0, U_plate_length/2.0, U_plate_thickness/2.0);
    G4LogicalVolume *Al_plate_L = new G4LogicalVolume(Al_plate_S, platemat, "Al_plate_L");

    U_plate_z = 1060.15; // det1z_actual + det1x_actual * tan(relative_angle);


    std::cout<<"\n\nPlacing Al plate at "<<U_plate_z<<"\n\n";
    // exit(0);
    G4PVPlacement *Al_plate_P = new G4PVPlacement(0,
                                    G4ThreeVector(0, 0, U_plate_z),
                                    Al_plate_L,
                                    "UraniumPlate",
                                    hall_L,
                                    false,
                                    0);

  }

  UVisAtt = new G4VisAttributes(G4Color(0,0,255));
  U_plate_L->SetVisAttributes(UVisAtt);

  // False near and far detectors right in front of the radiator and U plate to cover full space
  G4double F_plate_length = 18.0*inch;
  G4double F_plate_width  = 18.0*inch; // approximately
  G4double N_plate_length = 3*inch;
  G4double N_plate_width  = 3*inch; // approximately
  G4double F_plate_thickness = 0.1*mm; // thinnest DU plate we have is 1/32"

  G4Box *F_plate_S = new G4Box("F_plate_S", F_plate_width/2.0, F_plate_length/2.0, F_plate_thickness/2.0);
  G4Box *N_plate_S = new G4Box("N_plate_S", N_plate_width/2.0, N_plate_length/2.0, F_plate_thickness/2.0);

  G4LogicalVolume *F_plate_L = new G4LogicalVolume(F_plate_S, Air, "F_plate_L");
  G4LogicalVolume *N_plate_L = new G4LogicalVolume(N_plate_S, Air, "N_plate_L");

  //G4cout<<"\n\nFalse near: "<<collimator_lead_z-4.0*inch-1.5*mm<<"\n\n";

  if (incfalsedets)
  {
    G4PVPlacement *F_plate_P = new G4PVPlacement(0,
                                    G4ThreeVector(0, 0, U_plate_z-1.3*inch),
                                    F_plate_L,
                                    "FalseFar",
                                    hall_L,
                                    false,
                                    0);

    G4PVPlacement *N_plate_P = new G4PVPlacement(0,
                                    G4ThreeVector(0, 0, collimator_lead_z-4.0*inch-1.5*mm),
                                    N_plate_L,
                                    "FalseNear",
                                    hall_L,
                                    false,
                                    0);
  }

  //F_plate_L->SetSensitiveDetector(FarSD);
  //N_plate_L->SetSensitiveDetector(NearSD);

  BPVisAtt = new G4VisAttributes(G4Color(87./256, 193./256, 84./256));
  WVisAtt  = new G4VisAttributes(G4Color(200./256, 200./256, 200./256));


  // ------------------------------------------------------------------------------------
  // Multi-layer targets and foils
  // --maybe look into list initialization for this

  if (plateconf>3 && plateconf !=4)
  {
    // Thin genuine target: 09/13a, 09/14b, i.e. template I, i.e. NIM object 1
    G4Material* aMatSeptThinGenuine[] = { BoratedPoly, Al_mat,  DU_mat,  DU_mat,  DU_mat,  DU_mat,  Al_mat,  BoratedPoly };
    G4double    aThcSeptThinGenuine[] = { 3/4.0*inch,  0.25*mm, 0.91*mm, 0.95*mm, 0.94*mm, 0.92*mm, 0.25*mm, 3/4.0*inch  };
    vector<G4Material*> vMatSeptThinGenuine(std::begin(aMatSeptThinGenuine), std::end(aMatSeptThinGenuine));
    vector<G4double>    vThcSeptThinGenuine(std::begin(aThcSeptThinGenuine), std::end(aThcSeptThinGenuine));
    if (plateconf & 16384)
      BuildTarget("SeptThinGenuineTarget", vMatSeptThinGenuine, vThcSeptThinGenuine, 12.0*inch, 12.0*inch, "S");

    // Thin hoax target: 09/13b, 09/14a, i.e. hoax Ia, Ib, i.e. NIM object 2
    G4Material* aMatSeptThinHoax[] = { BoratedPoly, Al_mat,  Pb,          Pb,          Pb,          Pb,          Al_mat,  BoratedPoly };
    G4double    aThcSeptThinHoax[] = { 3/4.0*inch,  0.25*mm, 1/24.0*inch, 1/16.0*inch, 1/16.0*inch, 1/24.0*inch, 0.25*mm, 3/4.0*inch  };
    vector<G4Material*> vMatSeptThinHoax(std::begin(aMatSeptThinHoax), std::end(aMatSeptThinHoax));
    vector<G4double>    vThcSeptThinHoax(std::begin(aThcSeptThinHoax), std::end(aThcSeptThinHoax));
    if (plateconf & 8192)
      BuildTarget("SeptThinHoaxTarget", vMatSeptThinHoax, vThcSeptThinHoax, 12.0*inch, 12.0*inch, "S");

    // Thick genuine target: 09/15a, i.e. template II, i.e. NIM object 3
    G4Material* aMatSeptThickGenuine[] = { BoratedPoly, Al_mat,  DU_mat,  DU_mat,  DU_mat,  DU_mat,  DU_mat,  DU_mat,  DU_mat,  DU_mat,  Al_mat,  BoratedPoly };
    G4double    aThcSeptThickGenuine[] = { 3/4.0*inch,  0.25*mm, 0.91*mm, 0.95*mm, 0.94*mm, 0.92*mm, 0.89*mm, 0.85*mm, 0.88*mm, 0.85*mm, 0.25*mm, 3/4.0*inch  };
    vector<G4Material*> vMatSeptThickGenuine(std::begin(aMatSeptThickGenuine), std::end(aMatSeptThickGenuine));
    vector<G4double>    vThcSeptThickGenuine(std::begin(aThcSeptThickGenuine), std::end(aThcSeptThickGenuine));
    if (plateconf & 4096)
      BuildTarget("SeptThickGenuineTarget", vMatSeptThickGenuine, vThcSeptThickGenuine, 12.0*inch, 12.0*inch, "S");

    // Thick hoax target: 09/15d, i.e. hoax IIc, i.e. NIM object 4
    G4Material* aMatSeptThickHoax[] = { BoratedPoly, Al_mat,  Pb,          Pb,          Pb,          Pb,          Pb,          Pb,          Pb,          Pb,          Al_mat,  BoratedPoly };
    G4double    aThcSeptThickHoax[] = { 3/4.0*inch,  0.25*mm, 1/24.0*inch, 1/24.0*inch, 1/16.0*inch, 1/16.0*inch, 1/16.0*inch, 1/16.0*inch, 1/24.0*inch, 1/24.0*inch, 0.25*mm, 3/4.0*inch  };
    vector<G4Material*> vMatSeptThickHoax(std::begin(aMatSeptThickHoax), std::end(aMatSeptThickHoax));
    vector<G4double>    vThcSeptThickHoax(std::begin(aThcSeptThickHoax), std::end(aThcSeptThickHoax));
    if (plateconf & 2048)
      BuildTarget("SeptThickHoaxTarget", vMatSeptThickHoax, vThcSeptThickHoax, 12.0*inch, 12.0*inch, "S");

    // Partial hoax target: 09/15c, i.e. hoax IId, i.e. NIM object 5
    G4Material* aMatSeptPartialHoax[] = { BoratedPoly, Al_mat,  Pb,          Pb,          Pb,          Pb,          DU_mat,  DU_mat,  DU_mat,  DU_mat,  Al_mat,  BoratedPoly };
    G4double    aThcSeptPartialHoax[] = { 3/4.0*inch,  0.25*mm, 1/24.0*inch, 1/24.0*inch, 1/16.0*inch, 1/16.0*inch, 0.91*mm, 0.95*mm, 0.94*mm, 0.92*mm, 0.25*mm, 3/4.0*inch  };
    vector<G4Material*> vMatSeptPartialHoax(std::begin(aMatSeptPartialHoax), std::end(aMatSeptPartialHoax));
    vector<G4double>    vThcSeptPartialHoax(std::begin(aThcSeptPartialHoax), std::end(aThcSeptPartialHoax));
    if (plateconf & 1024)
      BuildTarget("SeptPartialHoaxTarget", vMatSeptPartialHoax, vThcSeptPartialHoax, 12.0*inch, 12.0*inch, "S");

    // standard 2.5" Al foil
    G4Material* aMatAlFoil[] = {Al_alloy};
    G4double    aThcAlFoil[] = {2.5*inch};
    vector<G4Material*> vMatAlFoil(std::begin(aMatAlFoil), std::end(aMatAlFoil));
    vector<G4double>    vThcAlFoil(std::begin(aThcAlFoil), std::end(aThcAlFoil));
    if (plateconf & 512)
    {
      U_plate_z = 1060.15;
      BuildTarget("AlFoil", vMatAlFoil, vThcAlFoil, 12.0*inch, 12.0*inch, "X");
    }

    // September nominal 4/32" DU foil, thicker in reality; add layers together so stackingAction works properly
    G4Material* aMatDUFoil[] = {DU_mat};
    G4double    aThcDUFoil[] = {(0.76+0.855+0.855+0.805)*mm};
    vector<G4Material*> vMatDUFoil(std::begin(aMatDUFoil), std::end(aMatDUFoil));
    vector<G4double>    vThcDUFoil(std::begin(aThcDUFoil), std::end(aThcDUFoil));
    if (plateconf & 256)
    {
      U_plate_z = 998.9;
      BuildTarget("DUFoil", vMatDUFoil, vThcDUFoil, 8.0*inch, 8.0*inch, "X");
    }

    // Load in Black Sea object for misalignment/rate tests on a spherical object
    if (plateconf & 128) {
      G4double r1_Pu_BS = 6.27*cm;
      G4double r2_Pu_BS = 6.70*cm;
      G4double r1_HE_BS = r2_Pu_BS;
      G4double r2_HE_BS = r1_HE_BS + 6.50*cm;
      G4double r1_WU_BS = r2_HE_BS;
      G4double r2_WU_BS = r1_WU_BS + 0.25*cm;

      G4Sphere *Pu_BS_S = new G4Sphere("Pu_BS_S", r1_Pu_BS, r2_Pu_BS, 0, twopi, 0, pi);
      G4Sphere *HE_BS_S = new G4Sphere("HE_BS_S", r1_HE_BS, r2_HE_BS, 0, twopi, 0, pi);
      G4Sphere *WU_BS_S = new G4Sphere("WU_BS_S", r1_WU_BS, r2_WU_BS, 0, twopi, 0, pi);

      G4LogicalVolume *Pu_BS_L = new G4LogicalVolume(Pu_BS_S, WgPu_simple, "Pu_BS_L"); // make sure to use WgPu_simple
      G4LogicalVolume *HE_BS_L = new G4LogicalVolume(HE_BS_S, HMX,         "HE_BS_L");
      G4LogicalVolume *WU_BS_L = new G4LogicalVolume(WU_BS_S, Uenriched95, "WU_BS_L");

      double BSz = collimator_lead_z+4*25.4+totalThickness/2.0+100.0*mm;

      G4VPhysicalVolume *Pu_BS_P = new G4PVPlacement(0, G4ThreeVector(0,0,BSz), Pu_BS_L, "BlackSeaPlutonium", hall_L, false, 0);
      G4VPhysicalVolume *HE_BS_P = new G4PVPlacement(0, G4ThreeVector(0,0,BSz), HE_BS_L, "BlackSeaHighExplosive", hall_L, false, 0);
      G4VPhysicalVolume *WU_BS_P = new G4PVPlacement(0, G4ThreeVector(0,0,BSz), WU_BS_L, "BlackSeaUranium", hall_L, false, 0);
    }
  }
}

  // ------------------------------------------------------------------------------------
  // List of physical volumes
  //
  double mTarget = 0.0;
  G4PhysicalVolumeStore *PVStore = G4PhysicalVolumeStore::GetInstance();
  //G4cout << G4endl;
  //G4cout << "List of physical volumes: " << G4endl;
  for (size_t i = 0; i < PVStore->size(); ++i)
  {
    G4VPhysicalVolume *pv = (*PVStore)[i];
    //G4cout << "  " << pv->GetName() << " (" << pv->GetLogicalVolume()->GetMaterial()->GetName() << ")" << G4endl;
    if (pv->GetName().contains("TargetLayer"))
      mTarget += pv->GetLogicalVolume()->GetMass();
  }
  //G4cout << G4endl;
  //G4cout << "Target mass = " << mTarget/kg << " kg." << G4endl;
  //G4cout << G4endl;

//  G4GDMLParser * parser = new G4GDMLParser();
//  parser->Write("ZKBrem1.gdml",hall_L);
//  exit(-9);

  return hall_P;
}


void geometryConstruction::BuildTarget(G4String name, vector<G4Material*> v_mats, vector<double> v_thickness, double width, double length, G4String opt)
{
  if (opt != "D" && opt != "X" && opt != "S") {//G4cout << "bad opt " << opt << ". Aborting..." << G4endl; exit(1);}

  const size_t nLayers = v_thickness.size();
  vector< G4Box* > v_boxes;
  vector< G4LogicalVolume* > v_log;

  // create the vector of solids and logical volumes
  double totalThickness = 0.0;
  for (size_t i = 0; i < nLayers; ++i){
    char buffS[128];
    char buffL[128];
    sprintf(buffS,"%s_%d_S",name.c_str(),(int)i);
    sprintf(buffL,"%s_%d_L",name.c_str(),(int)i);

    v_boxes.push_back((G4Box*) new G4Box(buffS, width/2.0, length/2.0, v_thickness[i]/2.0));

    v_log.push_back((G4LogicalVolume*) new G4LogicalVolume(v_boxes[i], v_mats[i], buffL));

    G4VisAttributes *layerVisAtt = new G4VisAttributes(); // default to grey
    if (v_mats[i] == U) layerVisAtt = UVisAtt;
    else if (v_mats[i] == BoratedPoly) layerVisAtt = BPVisAtt;
    v_log[i]->SetVisAttributes(layerVisAtt);

    totalThickness += v_thickness[i];
  }

  // create the mother solid and logical volumes
  char buffTS[128];
  char buffTL[128];
  sprintf(buffTS,"target%s_S",name.c_str());
  sprintf(buffTL,"target%s_L",name.c_str());
  G4Box *target_S = new G4Box(buffTS, width/2.0, length/2.0, totalThickness/2.0);
  G4LogicalVolume *target_L = new G4LogicalVolume(target_S, Vacuum, buffTL);
  target_L->SetVisAttributes(G4VisAttributes(false)); // invisible

  double target_x, target_y, target_z;
  target_x = 0*cm;
  target_y = 0*cm;
  target_z = (opt == "X" ? U_plate_z : 50.0*cm);
  target_z = (opt == "S" ? collimator_lead_z+4*25.4+totalThickness/2.0+20.0*mm : target_z);

  // place the mother volume
  G4VPhysicalVolume *target_P = new G4PVPlacement(0,
                                                  G4ThreeVector(target_x, target_y, target_z),
                                                  target_L,
                                                  name,
                                                  hall_L,
                                                  false,
                                                  0);

  // place the layers inside the mother volume
  double z_m = -totalThickness/2.0;
  for (size_t i = 0; i < nLayers; ++i){
    z_m += v_thickness[i]/2.0;
    char buffP[128];
    sprintf(buffP,"%sLayer%d",name.c_str(),(int)i);
    G4VPhysicalVolume *layer_P = new G4PVPlacement(0,
                                                   G4ThreeVector(0,0,z_m),
                                                   v_log[i],
                                                   buffP,
                                                   target_L,
                                                   false,
                                                   0);
    z_m += v_thickness[i]/2.0;
  }
}


G4LogicalVolume * geometryConstruction::BuildHPGe(unsigned int ind, G4ThreeVector pos, G4RotationMatrix * rot, double gA, double gB, double gC, double gD, double gE, double gF, double gG, double gH, double gI, double gJ, double gK, double gL, double gM, double gN)
{
  //G4cout<<"\n\nBuilding detector "<<ind<<"...\n";

  char buff[128];
  sprintf(buff,"%d",ind);
  G4String is(buff);

  // Crystal (subtracting radius of rounding from height)
  G4Tubs * HPGe_crystal_bulk_S = new G4Tubs(G4String("HPGe_crystal_bulk_S").append(is), 0.0, (gA-2.0*gM)/2.0, (gB - gJ)/2.0, 0.0, 2.0*M_PI);
  // Cap rounded corners
  G4Torus * HPGe_crystal_round_S = new G4Torus(G4String("HPGe_crystal_round_S").append(is),0.0, gJ-gM, (gA)/2.0-gJ, 0.0, 2.0*M_PI);
  // Cap torus fill in
  G4Tubs * HPGe_crystal_cap_S = new G4Tubs(G4String("HPGe_crystal_cap_S").append(is), 0.0, (gA-2.0*gM)/2.0-gJ-gM, (gJ-gM), 0.0, 2.0*M_PI);

  // Put together the bulk crystal shape (pre-bore)
  G4RotationMatrix * zerorot = new G4RotationMatrix(); zerorot->rotateX(0.0);
  G4UnionSolid * HPGe_crystal_c1_S = new G4UnionSolid(G4String("HPGe_crystal_c1_S").append(is),HPGe_crystal_bulk_S,HPGe_crystal_round_S, zerorot, G4ThreeVector(0.0,0.0,(gB - gJ)/2.0));
  G4UnionSolid * HPGe_crystal_c2_S = new G4UnionSolid(G4String("HPGe_crystal_c2_S").append(is),HPGe_crystal_c1_S,HPGe_crystal_cap_S, zerorot, G4ThreeVector(0.0,0.0,(gB - gJ)/2.0));
  G4UnionSolid * HPGe_crystal_c3_S = new G4UnionSolid(G4String("HPGe_crystal_c3_S").append(is),HPGe_crystal_round_S,HPGe_crystal_cap_S, zerorot, G4ThreeVector(0.0,0.0,0.0));

  // Put together the shape of the bore
  G4Tubs * bore_cyl_S = new G4Tubs(G4String("bore_cyl_S").append(is),0.0,gC/2.0,(gD-gC/2.0)/2.0,0.0,2.0*M_PI);
  G4Sphere * bore_cap_S = new G4Sphere(G4String("bore_cap_S").append(is),0.0,gC/2.0,0,2.0*M_PI,0,M_PI);
  G4UnionSolid * bore_S = new G4UnionSolid(G4String("bore_S").append(is),bore_cyl_S,bore_cap_S, zerorot, G4ThreeVector(0.0,0.0,(gD-gC/2.0)/2.0));

  // Make the bore
  G4SubtractionSolid * HPGe_crystal_S = new G4SubtractionSolid(G4String("HPGe_crystal_S").append(is),HPGe_crystal_c2_S,bore_S,zerorot,G4ThreeVector(0.0,0.0,-(gB - gJ)/2.0+(gD-gE)/2.0));

  // Make the main outer dead layer
  G4Tubs * HPGe_mDL_S = new G4Tubs(G4String("HPGe_mdL_S").append(is), (gA-2.0*gM)/2.0 , gA/2.0, (gB - gJ)/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume * hL = new G4LogicalVolume(HPGe_crystal_S, Ge, G4String("HPGe_crystal_L").append(is));
  // G4LogicalVolume * HPGe_mDL_L = new G4LogicalVolume(HPGe_mDL_S, Ge, G4String("HPGe_mDL_L").append(is));

  // Make the cap dead layer
  G4Torus * HPGe_cDL1_S = new G4Torus(G4String("HPGe_cDL1_S").append(is),0.0, gJ, (gA)/2.0-gJ, 0.0, 2.0*M_PI);
  G4Tubs * HPGe_cDL2_S = new G4Tubs(G4String("HPGe_cDL2_S").append(is), 0.0, (gA)/2.0-gJ, (gJ), 0.0, 2.0*M_PI);
  G4UnionSolid * HPGe_crystal_cDL3_S = new G4UnionSolid(G4String("HPGe_crystal_cDL3_S").append(is),HPGe_cDL1_S,HPGe_cDL2_S, zerorot, G4ThreeVector(0.0,0.0,0.0));
  G4SubtractionSolid * HPGe_crystal_cDL4_S = new G4SubtractionSolid(G4String("HPGe_crystal_cDL4_S").append(is),HPGe_crystal_cDL3_S,HPGe_crystal_c3_S, zerorot, G4ThreeVector(0.0,0.0,0.0));
  G4Box * HPGe_capcut_S = new G4Box(G4String("capcut_S").append(is), 500.0, 500.0, 20.0);
  G4SubtractionSolid * HPGe_crystal_capDL_S = new G4SubtractionSolid(G4String("HPGe_crystal_capDL_S").append(is),HPGe_crystal_cDL4_S,HPGe_capcut_S, zerorot, G4ThreeVector(0.0,0.0,-20.0));

  // Complete dead layer
  G4UnionSolid * HPGe_DL_S = new G4UnionSolid(G4String("HPGe_DL_S").append(is),HPGe_mDL_S,HPGe_crystal_capDL_S, zerorot, G4ThreeVector(0.0,0.0,(gB - gJ)/2.0));
  G4LogicalVolume * HPGe_DL_L = new G4LogicalVolume(HPGe_DL_S, Ge, G4String("HPGe_DL_L").append(is));

  // Inner can components
  G4Tubs * myl_cap_S = new G4Tubs(G4String("myl_cap_S").append(is),0.0,(90.00-2*gK)*mm/2.0,gH*mm/2.0,0.0,2.0*M_PI);
  G4LogicalVolume * myl_cap_L = new G4LogicalVolume(myl_cap_S, Mylar, G4String("myl_cap_L").append(is));

  G4Tubs * al_cap_S = new G4Tubs(G4String("al_cap_S").append(is),0.0,(90.00)*mm/2.0,gH*mm/2.0,0.0,2.0*M_PI);
  G4LogicalVolume * al_cap_L = new G4LogicalVolume(al_cap_S, Al, G4String("al_cap_L").append(is));

  G4Tubs * al_cani_S = new G4Tubs(G4String("al_cani_S").append(is),(90.00-2*gK)*mm/2.0,(90.00)*mm/2.0,(gF-3.0*mm-gH*mm)/2.0,0.0,2.0*M_PI);
  G4LogicalVolume * al_cani_L = new G4LogicalVolume(al_cani_S, Al, G4String("al_cani_L").append(is));

  G4Tubs * al_base_S = new G4Tubs(G4String("al_base_S").append(is),gC/2.0,(90.00)*mm/2.0,(3.0*mm)/2.0,0.0,2.0*M_PI);
  G4LogicalVolume * al_base_L = new G4LogicalVolume(al_base_S, Al, G4String("al_base_L").append(is));

  // Outer can
  G4Tubs * al_ocap_S = new G4Tubs(G4String("al_ocap_S").append(is),0.0,(108)*mm/2.0,gI*mm/2.0,0.0,2.0*M_PI);
  G4Tubs * al_ocan_S = new G4Tubs(G4String("al_ocan_S").append(is),(108-2*gI)*mm/2.0,(108)*mm/2.0,150*mm/2.0,0.0,2.0*M_PI);
  G4UnionSolid * al_ohouse_S = new G4UnionSolid(G4String("al_ohouse_S").append(is),al_ocan_S,al_ocap_S, 0, G4ThreeVector(0,0,150/2.0));
  G4LogicalVolume * al_ohouse_L = new G4LogicalVolume(al_ohouse_S, Al, G4String("al_ohouse_L").append(is));

  // Containment boxes for each detector
  G4Box * det_box_S = new G4Box(G4String("det_box_S").append(is), 108.001/2.0, 108.001/2.0, 230.0/2.0);
  G4LogicalVolume *det_box1_L = new G4LogicalVolume(det_box_S, Air, G4String("det_box_L").append(is));

  G4PVPlacement * det_box1_P = new G4PVPlacement(rot,
                                                    pos,
                                                    det_box1_L,
                                                    G4String("detbox").append(is),
                                                    hall_L,
                                                    false,
                                                    0);

//G4ThreeVector(det1x_actual-6.3*sin(relative_angle)-5.0*cos(relative_angle) + 3.0*inch*cos(relative_angle), det1y_actual-6.3, det1z_actual
//                                                             - 6.3*cos(relative_angle) + 5.0*sin(relative_angle) -3.0*inch*sin(relative_angle))

  // HPGe detector closest to radiator, bottom row, with its can
  G4PVPlacement *HPGe_crystal_1_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,0+59.0),
                                                    hL,
                                                    G4String("HPGeCrystal").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);

  G4PVPlacement *HPGe_DL_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,0+59.0),
                                                    HPGe_DL_L,
                                                    G4String("HPGe_DL_").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);

  G4PVPlacement *myl_cap1_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,gB/2+gH/2+gJ/2.0+59.0),
                                                    myl_cap_L,
                                                    G4String("MylarCap").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);
  G4PVPlacement *al_cap1_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,gB/2+gH/2+gH+gJ/2.0+59.0),
                                                    al_cap_L,
                                                    G4String("AlCap").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);
  G4PVPlacement *al_cani1_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,gB/2+gH-(gF-3.0*mm-gH*mm)/2+gJ/2.0+59.0),
                                                    al_cani_L,
                                                    G4String("AlCanInner").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);

  G4PVPlacement *al_base1_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,gB/2+gH-(gF-3.0*mm-gH*mm)/2+gJ/2-(gF-3.0*mm-gH*mm)/2-3.0/2.0+59.0),
                                                    al_base_L,
                                                    G4String("AlBase").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);
  G4PVPlacement *al_ohouse_P = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,-30+59.0+(gB-69.4)/2),
                                                    al_ohouse_L,
                                                    G4String("AlHouse").append(is),
                                                    det_box1_L,
                                                    false,
                                                    0);

  //G4cout<<"Detector "<<ind<<" complete.\n\n";

  G4ThreeVector foilToDetRay = pos + HPGe_crystal_1_P->GetTranslation() - G4ThreeVector(0,0,U_plate_z);
  G4ThreeVector beamDir = G4ThreeVector(0,0,1);
  G4double thetaAngle = foilToDetRay.angle(beamDir)*180.0/3.1415;
  //G4cout << "thetaAngle = " << thetaAngle << " degrees." << G4endl;

  return hL;
}

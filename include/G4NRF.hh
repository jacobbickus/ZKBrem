//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NRF.hh,v 1.1.1.1 2007/01/13 19:44:02 jordan Exp $
// GEANT4 tag $Name:  $
//
//      ------------ G4NRF physics process --------
//
//      File name: G4NRF
//      Author:        David Jordan (david.jordan@pnl.gov)
//      Creation date: October 2006
//
//      How to use this physics process:
//
//      (1) Store ENSDF-derived nuclear level and gamma emission
//          files in a directory pointed to by the environment 
//          variable $G4NRFGAMMADATA.  These files include:
//             z[Z-value].a[A-value] -- default photon evaporation files
//             gamma_table_nnn.dat   -- supplementary info for gamma emission
//             level_table_nnn.dat   -- supplementary level info
//             ground_state_properties.dat
//
//      (2) Place the following code in the PhysicsList class:
//
//       #include "G4NRF.hh"
//       theParticleIterator->reset();
//       while( (*theParticleIterator)() ){
//         G4ParticleDefinition* particle = theParticleIterator->value();
//         G4ProcessManager* pmanager     = particle->GetProcessManager();
//         G4String particleName          = particle->GetParticleName();
//
//         if (particleName == "gamma") {
//             .... other gamma processes ...
//
//           pmanager->AddDiscreteProcess(new G4NRF("NRF"));
//
//             .... etc. (other gamma processes) ...
//             
//          }
//       }
//
//       To turn on verbose output from NRF package (diagnostic info):
//       
//         G4bool Verbose;
//         pmanager->AddDiscreteProcess(new G4NRF("NRF", Verbose=true));
//
//       Alternatively, declare a pointer to the NRF process and call the
//       SetVerbose() member function:
//
//       G4NRF* pNRF = new G4NRF("NRF");
//       pmanager->AddDiscreteProcess(pNRF);
//       pNRF->SetVerbose();
//
//      (3) Optional features:
//
//       To turn on verbose output from NRF package (diagnostic info):
//       
//         G4bool Verbose;
//         pmanager->AddDiscreteProcess(new G4NRF("NRF", Verbose=true));
//
//       Alternatively, declare a pointer to the NRF process and call the
//       SetVerbose() member function:
//
//       G4NRF* pNRF = new G4NRF("NRF");
//       pmanager->AddDiscreteProcess(pNRF);
//       pNRF->SetVerbose();
//
//       To disable gamma-gamma angular correlation: 
//
//       pNRF->ForceIsotropicAngCor();
//
//      (4) Material definition requirements: Materials must be constructed
//          from elements containing G4Isotopes in order to trigger 
//          calculation of NRF cross sections.  (Otherwise no NRF interaction
//          will occur.)
//
//          Example (in DetectorConstruction class):
//
//        // Natural Oxygen
//
//        G4Isotope *O16 = new G4Isotope(name="O16", iz=8, n=16,  a=15.995*g/mole);
//        G4Isotope *O17 = new G4Isotope(name="O17", iz=8, n=17,  a=16.999*g/mole);
//        G4Isotope *O18 = new G4Isotope(name="O18", iz=8, n=18, a=17.992*g/mole);
//        G4Element *NatO = new G4Element
//                     (name="Natural Oxygen", symbol="O", ncomponents=3);
//        NatO->AddIsotope(O16, abundance=99.757*perCent);
//        NatO->AddIsotope(O17, abundance= 0.038*perCent);
//        NatO->AddIsotope(O18, abundance= 0.205*perCent);
// 
// ------------------------------------------------------------

// G4NRF Class description
//
// This class manages NRF
// it inherites from G4VDiscreteProcess
//
// Class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4NRF_h
#define G4NRF_h 1

#include "G4ios.hh" 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OrderedTable.hh" 
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4NRFNuclearLevelManager.hh"
#include "AngularCorrelation.hh"

#include "G4VProcess.hh"

#include "G4Integrator.hh"

#include <fstream>
using std::ofstream;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
class G4NRF : public G4VDiscreteProcess
 
{ 
  public:
 
     G4NRF(const G4String& processName = "NRF", G4bool Verbose_in = false);
 
    ~G4NRF();

     G4bool IsApplicable(const G4ParticleDefinition&);
     
     void PrintInfoDefinition();
           
     G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition );
 
     G4VParticleChange *PostStepDoIt(const G4Track& track,         
                                     const G4Step&  step);                 
     void SetVerbose() {Verbose = true; };
     void ForceIsotropicAngCor() {ForceIsotropicAngularCorrelation = true; };

     void SetParam_x(G4double x);
     G4double GetParam_x() const;

     void SetParam_t(G4double t);
     G4double GetParam_t() const;

     inline G4double fast_exp(G4double y) const;

     G4double expIntegrand(G4double y) const;
     G4double HermiteIntegrand(G4double a) const;

     G4double PsiIntegral(G4double x, G4double t);

     void BuildNumIntTable();
     std::vector< std::vector<G4double> > num_int_table;
     G4double PsiIntegralLookup(G4double x, G4double t);

     void print_to_standalone(ofstream& file);

  private:
  
     G4NRF & operator=(const G4NRF &right);
     G4NRF(const G4NRF&);
     G4NRFNuclearLevelManager* pNuclearLevelManager;
     G4double NRF_xsec_calc_legacy(G4double GammaEnergy, G4double J0, G4int A, const G4NRFNuclearLevel* pLevel);
     G4double NRF_xsec_calc(G4double GammaEnergy, G4double J0, G4int A, const G4NRFNuclearLevel* pLevel);

     void SetupMultipolarityInfo(const G4int nLevel, 
				 const G4double E_gamma,
				 const G4int jgamma, 
				 const G4NRFNuclearLevel* pLevel,
				 const G4NRFNuclearLevel* pLevel_next,
				 G4double& J0, G4double& J, G4double& Jf,
				 G4int& L1, G4int& L2,
				 G4double& Delta1, G4double& Delta2);

    void AssignMultipoles(const G4NRFNuclearLevel* pLevel, 
			  const G4double E_gamma,
			  const G4int jgamma,
			  const G4double Ji, const G4double Pi,
			  const G4double Jf, const G4double Pf,
			  G4int& L, G4double& Delta);

     G4ThreeVector SampleCorrelation(const G4double Ji, const G4double J,  const G4double Jf, 
				     const G4int L1, const G4int L2, 
				     const G4double Delta1, const G4double Delta2);

     G4ThreeVector SampleIsotropic();

     G4int FindMin_L(const G4double Ji, const G4double Pi, 
		     const G4double Jf, const G4double Pf, char& transition);

     G4bool ForbiddenTransition(const G4NRFNuclearLevel* pLevel,
				const G4NRFNuclearLevel* pLevel_next);
  
     G4double MixingRatio_WeisskopfEstimate(char multipole, const G4int L,
					    const G4double E_gamma,
					    const G4double Pi,
					    const G4double Pf);

     G4double Lamda_Weisskopf(char multipole, const G4int L, 
			      const G4double E_gamma);

     G4int A_excited;
     G4int Z_excited;
     const G4NRFNuclearLevel* pLevel_excited;

     Angular_Correlation* pAngular_Correlation;

     G4bool Verbose;
     G4bool ForceIsotropicAngularCorrelation;

     G4double param_x;
     G4double param_t;

     // numerical integration
     G4Integrator<const G4NRF, G4double(G4NRF::*)(G4double) const> integrator;
     G4int nIterations;
     G4double ztol2;

     // declare parameters for numerical integration table
     G4int nt;
     G4double tmin;
     G4double tmax;
     G4double dt;

     G4int nx;
     G4double xmin;
     G4double xmax;
     G4double dx;
};

inline G4bool G4NRF::IsApplicable(
                            const G4ParticleDefinition& particle)
{
  return (&particle == G4Gamma::Gamma());
}

#endif


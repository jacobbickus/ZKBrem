#ifndef G4NRFPhysics_h
#define G4NRFPhysics_h 1

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsListHelper.hh"
#include "G4LossTableManager.hh"
#include "G4BuilderType.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

#include "G4Gamma.hh"

#include "G4NRF.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4VModularPhysicsList.hh"


class G4NRFPhysics : public G4VPhysicsConstructor
{
	public:
		G4NRFPhysics(const G4String &name = "NRF");

		virtual ~G4NRFPhysics();

		virtual void ConstructParticle();
    virtual void ConstructProcess();

};

#endif
#ifndef particleSplittingOperation_hh
#define particleSplittingOperation_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForNothing.hh"

class G4LogicalVolume;

class particleSplittingOperation : public G4VBiasingOperation
{
public:
	particleSplittingOperation(G4String name);
	virtual ~particleSplittingOperation();

	virtual G4double DistanceToApplyOperation(const G4Track*, G4double, G4ForceCondition* condition);
	virtual G4VParticleChange *GenerateBiasingFinalState(const G4Track*, const G4Step*);

	inline void SetSplittingFactor(G4int splittingFactor) {fSplittingFactor = splittingFactor;};
	inline G4int GetSplittingFactor() {return fSplittingFactor;};

	// unused but required virtual methods
	virtual const G4VBiasingInteractionLaw *ProvideOccurenceBiasingInteractionLaw(const G4BiasingProcessInterface* /*, G4ForceCondition&*/) {return 0;};
	virtual G4VParticleChange *ApplyFinalStateBiasing(const G4BiasingProcessInterface*, const G4Track*, const G4Step* /*, G4bool&*/) {return 0;};

private:
	G4ParticleChange fParticleChange;
	G4ParticleChangeForNothing fParticleChangeForNothing;
	G4int fSplittingFactor;
};

#endif

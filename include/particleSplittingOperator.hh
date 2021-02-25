#ifndef particleSplittingOperator_hh
#define particleSplittingOperator_hh 1

#include "G4VBiasingOperator.hh"
#include "particleSplittingOperation.hh"

class particleSplittingOperator : public G4VBiasingOperator
{
public:
	particleSplittingOperator();
	virtual ~particleSplittingOperator();

	particleSplittingOperation *GetSplitAndKillOperation() const {return fSplitAndKillOperation;};


private:
	particleSplittingOperation *fSplitAndKillOperation;
	G4int fSplittingFactor;

	// used for splitting/killing
	virtual G4VBiasingOperation *ProposeNonPhysicsBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess);

	// unused, but required virtual functions
	virtual G4VBiasingOperation *ProposeOccurenceBiasingOperation(const G4Track*, const G4BiasingProcessInterface* ) { return 0; };
  virtual G4VBiasingOperation *ProposeFinalStateBiasingOperation(const G4Track*, const G4BiasingProcessInterface* ) { return 0; };

};

#endif

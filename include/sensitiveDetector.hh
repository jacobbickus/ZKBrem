#ifndef sensitiveDetector_h
#define sensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include <vector>

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;
class rootStorageManager;

class sensitiveDetector : public G4VSensitiveDetector
{
  public:
    sensitiveDetector(const G4String&);
    virtual ~sensitiveDetector();

    void   Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    sensitiveDetector & operator=(const sensitiveDetector &right);
    sensitiveDetector(const sensitiveDetector&);
    
    G4double totalEne;
    G4double trackWeight;
    G4double gE, gt, gx, gy;
    G4bool qNRF1, qNRF2, qNRF3;

    bool pdone;

    rootStorageManager *rm;
    
};

#endif

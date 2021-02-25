#ifndef eventInformation_h
#define eventInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4Allocator.hh"
#include "G4VUserEventInformation.hh"
#include "TObject.h"

class eventInformation : public G4VUserEventInformation
{
	public:
		eventInformation();
		eventInformation(const G4Event*);
		eventInformation(const eventInformation*);
		virtual ~eventInformation();

		inline G4double GetWeight() const {return weight;}
		void SetWeight(G4double);

		inline G4double GetBeamEnergy() const {return beamEnergy;}
		void SetBeamEnergy(G4double);
		
		inline G4double GetBeamAngle() const {return beamAngle;}
		void SetBeamAngle(G4double);

		inline G4double GetGammaEnergy() const {return gammaEnergy;}
		void SetGammaEnergy(G4double);

		inline G4double GetVertexX() const {return vertexX;}
		void SetVertexX(G4double);

		inline G4double GetVertexY() const {return vertexY;}
		void SetVertexY(G4double);

		inline G4double GetVertexZ() const {return vertexZ;}
		void SetVertexZ(G4double);

		inline G4double GetVertexTheta() const {return vertexTheta;}
		void SetVertexTheta(G4double);

		inline G4double GetVertexPhi() const {return vertexPhi;}
		void SetVertexPhi(G4double);

		inline G4double GetPlateDepartureE() const {return platedE;}
		void SetPlateDepartureE(G4double);

		inline G4ThreeVector GetFoilTrackInit() const {return foilTrackInit;}
		void SetFoilTrackInit(G4ThreeVector);

		inline G4double GetFoilPathLength() const {return foilPathLength;}
		void SetFoilPathLength(G4double);

		inline ULong64_t GetGeneratorIndex() const {return genind;}
		void SetGeneratorIndex(ULong64_t x);

		void Print() const;

		G4double gEnear;
		G4double gEfar;

	private:
		G4double weight;
		G4double beamEnergy;
		G4double gammaEnergy;
		G4double beamAngle;
    G4double vertexX;
    G4double vertexY;
    G4double vertexZ;
    G4double vertexTheta;
    G4double vertexPhi;
    G4double platedE;
    G4ThreeVector foilTrackInit;
    G4double foilPathLength;
    ULong64_t genind;
};

#endif

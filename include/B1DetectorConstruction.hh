#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4LogicalVolume* GetFrontChVolume() const {return fFrontChVolume;}

  protected:
    G4LogicalVolume*  fScoringVolume;
    G4LogicalVolume* fFrontChVolume;
    std::vector<G4LogicalVolume*> vLayers;
};


#endif

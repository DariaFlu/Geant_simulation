#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"


B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }


B1DetectorConstruction::~B1DetectorConstruction()
{ }


G4VPhysicalVolume* B1DetectorConstruction::Construct()


{

    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();

    // Envelope parameters
    G4double env_sizeX = 100*mm, env_sizeY = 100*mm, env_sizeZ = 2000*mm;
    G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = true;

    //Create the box - mother volume "world"
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4double world_x = 100.*mm;
    G4double world_y = 100.*mm;
    G4double world_z = 2000.*mm;

    G4Box* solidWorld =
        new G4Box("world_box",
                    0.5*world_x, 0.5*world_y, 0.5*world_z);

    G4LogicalVolume* logicWorld =
        new G4LogicalVolume(solidWorld,
                            world_mat,
                            "World");

    G4VPhysicalVolume* physWorld =
        new G4PVPlacement(0,
                          G4ThreeVector(),
                          logicWorld,
                          "World",
                          0,
                          false,
                          0,
                          checkOverlaps);

    // Envelope
    G4Box* solidEnv =
        new G4Box("Envelope",
                   0.5*env_sizeX, 0.5*env_sizeY, 0.5*env_sizeZ);

    G4LogicalVolume* logicEnv =
        new G4LogicalVolume(solidEnv,
                            env_mat,
                            "Envelope");

        new G4PVPlacement(0,
                          G4ThreeVector(),
                          logicEnv,
                          "Envelope",
                          logicWorld,
                          false,
                          0,
                          checkOverlaps);

/*
The beam pipe exit window is 55 micrometers  thick
the beam pipe exit window is from havar alloy   8.3 g/cm3
  Havar is composed of 42.0% (41-44%) of cobalt,
                     19.5% (19-21%) of chromium,
                     12.7% (12-14%) of nickel,
                     2.7% (2.3-3.3%) of tungsten,
                     2.2% (2-2.8%) of molybdenum,
                     1.6% (1.35-1.8%) of manganese,
                     0.2% (0.17-0.23%) of carbon,
                     0.02-0.08% of beryllium, and balance of iron.
*/

    G4double densityHavar = 8300.* kg/m3;
    std::vector<G4String> Havar_elm;
    std::vector<G4double> Havar_weight;
    G4double temperature = 300.; //K->0C 26.85 celsius
    G4double pressure    = 10^5; //Pa

    Havar_elm.push_back("Co");
    Havar_elm.push_back("Cr");
    Havar_elm.push_back("Ni");
    Havar_elm.push_back("W");
    Havar_elm.push_back("Mo");
    Havar_elm.push_back("Mn");
    Havar_elm.push_back("C");
    Havar_elm.push_back("Be");
    Havar_elm.push_back("Fe");

    Havar_weight.push_back(0.420);
    Havar_weight.push_back(0.195);
    Havar_weight.push_back(0.127);
    Havar_weight.push_back(0.027);
    Havar_weight.push_back(0.022);
    Havar_weight.push_back(0.016);
    Havar_weight.push_back(0.002);
    Havar_weight.push_back(0.0005);
    Havar_weight.push_back(0.1905);

    G4Material * HavarAlloy = nist->ConstructNewMaterial ("HavarAlloy",
                                                           Havar_elm,
                                                           Havar_weight,
                                                           densityHavar,
                                                           true,
                                                           kStateSolid,
                                                           temperature,
                                                           pressure
                                                         );

    G4Material* HavarAlloy_mat = nist->FindOrBuildMaterial("HavarAlloy");

    G4double thick  = 55.*um;
    G4double radius = 1.5*mm;

    /*G4Box* beamWindow =
        new G4Box("beamWindow",
                  thick,
                  thick,
                  thick
        );*/
    G4Tubs* beamWindow =
        new G4Tubs ( "beamWindow", 0, radius, thick, 0.*deg, 360.*deg);

    G4LogicalVolume* beamWindow_log =
        new G4LogicalVolume(beamWindow, HavarAlloy_mat,"beamWindow_log"); //The beam pipe exit window

    //Create a Physical Volume
    G4double windowPos_x = 0.0*mm;
    G4double windowPos_y = 0.0*mm;
    G4double windowPos_z = -400.0*mm;

    new G4PVPlacement(0 , // no rotation
        G4ThreeVector(windowPos_x, windowPos_y, windowPos_z),
                    beamWindow_log,
                    "beamWindowPlace",
                    logicWorld, //LOGICAL VOLUME
                    false,
                    0);

    //Create the ionization chamber
    G4double innerRadiusOfTheChamber   = 0.*mm;
    G4double outerRadiusOfTheChamber   = 6.95*mm;
    G4double hightOfTheChamber         = 23.6*mm;
    //The wall of the graphite cylinder has a thickness of 0.09 mm.
    G4double innerRadiusOfTheGraphite  = 6.1*mm;
    G4double outerRadiusOfTheGraphite  = 6.28*mm;
    //The wall of the PMMA cylinder is 0.335mm thick.
    G4double innerRadiusOfThePMMA      = 6.28*mm;
    //The diameter of the central aluminum electrode is 1.1 mm
    G4double outerRadiusOfTheRod       = 1.1*mm;
    G4double hightOfTheRod             = 21.2*mm;
    //AIR
    G4double hightOfAir                = 23.0*mm;
    //Angle
    G4double startAngleOfTheChamber    = 0.*deg;
    G4double spanningAngleOfTheChamber = 360.*deg;

    G4Tubs* alRod =
        new G4Tubs ("alRod",
            innerRadiusOfTheChamber*0.5,
            outerRadiusOfTheRod*0.5,
            hightOfTheRod*0.5,
            startAngleOfTheChamber,
            spanningAngleOfTheChamber );
    G4Material* alRod_mat = nist->FindOrBuildMaterial("G4_Al");

    G4Tubs* air =
        new G4Tubs ("air",
            innerRadiusOfTheChamber*0.5,
            innerRadiusOfTheGraphite*0.5,
            hightOfAir*0.5,
            startAngleOfTheChamber,
            spanningAngleOfTheChamber
        );
    G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4Tubs* graphite =
        new G4Tubs ("graphite",
            innerRadiusOfTheChamber*0.5,
            outerRadiusOfTheGraphite*0.5,
            hightOfAir*0.5,
            startAngleOfTheChamber,
            spanningAngleOfTheChamber
        );
    G4Material* graphite_mat = nist->FindOrBuildMaterial("G4_GRAPHITE");

    G4Tubs* pmma =
        new G4Tubs ("pmma",
            innerRadiusOfTheChamber*0.5,
            outerRadiusOfTheChamber*0.5,
            hightOfTheChamber*0.5,
            startAngleOfTheChamber,
            spanningAngleOfTheChamber
        );

    G4double density;
    G4int polyPMMA = 1;

    std::vector<G4String> PMMA_elm;
    std::vector<G4int> PMMA_nbAtoms;

    PMMA_elm.push_back("H");
    PMMA_elm.push_back("C");
    PMMA_elm.push_back("O");
    PMMA_nbAtoms.push_back(3+2*polyPMMA);
    PMMA_nbAtoms.push_back(6+2*polyPMMA);
    PMMA_nbAtoms.push_back(2);

    G4Material* PMMA     = nist->ConstructNewMaterial("PMMA",PMMA_elm,PMMA_nbAtoms, density=1190*kg/m3);
    G4Material* pmma_mat = nist->FindOrBuildMaterial("PMMA");

    /*G4Tubs* tracker_Chamber =
        new G4Tubs("tracker_tube",
                    innerRadiusOfTheChamber,
                    outerRadiusOfTheChamber,
                    hightOfTheChamber,
                    startAngleOfTheChamber,
                    spanningAngleOfTheChamber);*/


    //Logical volume

    G4LogicalVolume* alRod_log =
        new G4LogicalVolume(alRod, alRod_mat,"alRod_log"); //Aluminum electrode

    G4LogicalVolume* air_log =
        new G4LogicalVolume(air, air_mat,"air_log"); //Air

    G4LogicalVolume* graphite_log =
        new G4LogicalVolume(graphite, graphite_mat,"graphite_log"); //graphite

    G4LogicalVolume* pmma_log =
        new G4LogicalVolume(pmma, pmma_mat, "pmma_log"); //


    //Create a Physical Volume
    G4double chamberPos_x = 0.0*mm;
    G4double chamberPos_y = 0.0*mm;
    G4double chamberPos_z = 1800.0*mm;

    G4RotationMatrix *rotat = new G4RotationMatrix;
    rotat->rotateX(90*deg);
    rotat->rotateY(180*deg);
    rotat->rotateZ(0*deg);



    new G4PVPlacement(rotat, // rotation
        G4ThreeVector(chamberPos_x*0.5, chamberPos_y*0.5, chamberPos_z*0.5),
                    pmma_log,
                    "pmmPlace",
                    logicWorld, //LOGICAL VOLUME
                    false,
                    0);

    new G4PVPlacement(rotat, // rotation
        G4ThreeVector(chamberPos_x*0.5, chamberPos_y*0.5, chamberPos_z*0.5),
                    graphite_log,
                    "graphitePlace",
                    logicWorld, //LOGICAL VOLUME
                    false,
                    0);

    new G4PVPlacement(rotat, //rotation
        G4ThreeVector(chamberPos_x*0.5, chamberPos_y*0.5, chamberPos_z*0.5),
                    air_log,
                    "airPlace",
                    logicWorld, //LOGICAL VOLUME
                    false,
                    0);

    new G4PVPlacement(rotat, // no rotation
        G4ThreeVector(chamberPos_x*0.5, chamberPos_y*0.5, chamberPos_z*0.5),
                    alRod_log,
                    "alRodPlace",
                    logicWorld, //LOGICAL VOLUME
                    false,
                    0);

    fScoringVolume = air_log;
    return physWorld;
}

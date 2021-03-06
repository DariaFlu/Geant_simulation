#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <vector>


B1RunAction::B1RunAction() : G4UserRunAction(),
  fEdep(0.),  //Step is used for ionization chamber
  fEdep2(0.), //Step2 is used for ionization chamber
  fEbef(0.), //The energy before proton hitting chamber
  fNofLayers(100), //Number of layers
  hist(500,0), //histogram
  ftotEnerg(0.)
{
  // add new units for dose
  //
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
  vDepoEnr_run.resize(fNofLayers, 0.);
}


B1RunAction::~B1RunAction()
{}


void B1RunAction::BeginOfRunAction(const G4Run*)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  ftotEnerg = 0.;
  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  for(size_t i = 0; i<500; i++){ hist[i]=0; }
    vDepoEnr_run.resize(fNofLayers, 0.);
  fEbef = 0.;
  nHitFCBef = 0;
}


void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
        = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : "
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;


     G4cout << "Writing spectrum in spectrum.dat file..." << G4endl;
     std::ofstream filed0("spectrum.dat", std::ios::trunc);
     double bin_width = (HIST_MAX - HIST_MIN) / NOBINS;
     for(int i = 0; i<NOBINS; i++){
         double energy0 = i*bin_width + HIST_MIN;
         filed0 << energy0 << " " << hist[i] << std::endl;
     }
     filed0.close();

     G4cout << "Writing spectrum in layers.dat file..." << G4endl;
     std::ofstream filed1("layers.dat", std::ios::trunc);
     //double bin_width = (HIST_MAX - HIST_MIN) / NOBINS;
     for(int i = 0; i<fNofLayers; i++){
         //double energy0 = i*bin_width + HIST_MIN;
         filed1 << i << " " << vDepoEnr_run[i]/MeV/nofEvents << std::endl;
     }
     filed1.close();

     G4cout<<"Energy before chamber: " <<fEbef/MeV/nHitFCBef<<G4endl;
     G4cout<<"Total energy in chamber: " <<ftotEnerg/MeV/nHitFCBef<<G4endl;
     G4cout << "Done!" << G4endl;
 }

void B1RunAction::AddEdep(G4double edep){
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void B1RunAction::PutInHisto(G4double edep){
    double bin_width = (HIST_MAX - HIST_MIN) / NOBINS;
    int index0 = int(floor((edep-HIST_MIN)/bin_width/keV)); //energy spread simulation is used
    if(index0 > 0 && index0 < NOBINS){ hist[index0]++; }
    ftotEnerg+=edep;
}

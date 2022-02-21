#pragma once

#ifndef B1SteppingAction_h
#define B1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "B1RunAction.hh"

#include <string>
#include <vector>

class B1EventAction;
class B1RunAction;

class G4LogicalVolume;

/// Stepping action class
///

class B1SteppingAction : public G4UserSteppingAction
{
  public:
    B1SteppingAction(B1EventAction* eventAction, B1RunAction* runAction);
    virtual ~B1SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    B1EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
    B1RunAction*    fRunAction;
};

//CHECKS IF VOLUME NAME IS A NUMBER
bool is_number(const std::string& stn);


#endif

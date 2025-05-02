#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "detector.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include <vector>
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Scintillation.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    virtual G4VPhysicalVolume * Construct();

    G4LogicalVolume *GetScoringVolume() const {return fScoringVolume;} 

private:
    G4LogicalVolume * logicDetector;
    virtual void ConstructSDandField();
    
    // Define all the variables we need
    G4int nCols, nRows;
    G4Box *solidWorld; 
    G4Tubs *solidXenon, *solidPhotocathode, *solidDetector, *solidPMT;
    G4LogicalVolume *logicWorld, *logicPhotocathode, *logicXenon, *logicPMT;
    G4VPhysicalVolume *physWorld, *physDetector, *physXenon, *physPMT_Top, *physPMT_Bottom;
    G4String name, symbol;
    G4double z, a, density, temperature, pressure, pmtRadius;
    G4int ncomponents;
    G4Material *lXeMat, *gXeMat, *worldMat, *quartzMat, *photocathodeMat, *PTFEMat, *cryoMat, *copperMat, *kovarMat;
    std::vector<G4ThreeVector> pmtPositions;

    void DefineMaterial();

    void DefinePMTs();

    G4GenericMessenger *fMessenger;

    G4LogicalVolume *fScoringVolume;

    // Implement quantum efficiency for measured PMT spectrum
    G4MaterialPropertiesTable* LoadQuantumEfficiency(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            G4cerr << "Error: Could not open file " << filename << G4endl;
            return nullptr;
        }

        std::vector<G4double> energies;
        std::vector<G4double> efficiencies;

        G4double energy, efficiency;
        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            if (iss >> energy >> efficiency) {
                energies.push_back(energy);
                efficiencies.push_back(efficiency);
            }
        }
        
        file.close();

        if (energies.empty()) {
            G4cerr << "Error: No data found in " << filename << G4endl;
            return nullptr;
        }

        // Create Material Properties Table for Photocathode
        G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("EFFICIENCY", energies.data(), efficiencies.data(), energies.size());
        std::vector<G4double> cathmetal_PP   = { 5.0*eV, 6.69*eV, 7.50*eV };
        std::vector<G4double> cathmetal_RIND = { 1.51, 1.57, 1.61 };     // ref index
        std::vector<G4double> cathmetal_ABSL = { 1.e-20*m,  1.e-20*m,  1.e-20*m };// atten length
        mpt->AddProperty("REFLECTIVITY", cathmetal_PP, {0,0,0});
        mpt->AddProperty("RINDEX", cathmetal_PP, cathmetal_RIND);
        mpt->AddProperty("ABSLENGTH", cathmetal_PP, cathmetal_ABSL);

        return mpt;
    }


    
};

#endif

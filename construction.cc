#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{   

    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");


    DefineMaterial();

}

void MyDetectorConstruction::DefineMaterial()
{
    G4NistManager *nist = G4NistManager::Instance();

   /* G4Material* LXe = new G4Material("LXe", 54, 131.29 * g/mole,
                                      3.02 * g/cm3, kStateLiquid,
                                      173.15 * kelvin, 1.5 * atmosphere);

    std::vector<G4double> lxe_Energy   = { 7.0*eV , 7.07*eV, 7.14*eV };
    std::vector<G4double> lxe_SCINT = { 0.1, 1.0, 0.1 };
    std::vector<G4double> lxe_RIND  = { 1.59 , 1.57, 1.54 };
    std::vector<G4double> lxe_ABSL  = { 35.*cm, 35.*cm, 35.*cm}; //atten length

    G4MaterialPropertiesTable* fLXe_mt = new G4MaterialPropertiesTable();
    fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", lxe_Energy, lxe_SCINT);
    fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", lxe_Energy, lxe_SCINT);
    fLXe_mt->AddProperty("RINDEX", lxe_Energy, lxe_RIND);
    fLXe_mt->AddProperty("ABSLENGTH", lxe_Energy, lxe_ABSL);
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000./MeV);
    fLXe_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20.*ns);
    fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45.*ns);
    LXe->SetMaterialPropertiesTable(fLXe_mt);
    // Set the Birks Constant for the LXe scintillator
    LXe->GetIonisation()->SetBirksConstant(0.126*mm/MeV);   */
    
    
    G4Element* elementXe = new G4Element( "Xenon", "Xe", 54., 131.29*g/mole );
    G4Material* LXe = new G4Material
        ("LXe", 3.02*g/cm3, 1, kStateLiquid, 173.15*kelvin, 1.5*atmosphere );
    G4Material* GXe = new G4Material
        ("GXe", 0.005887*g/cm3, 1, kStateGas, 173.15*kelvin, 1.5*atmosphere );
    LXe->AddElement( elementXe, 1);
    GXe->AddElement( elementXe, 1);

    //  std::vector<G4double> LXe_PP    = { 7.07*eV, 7.07*eV };
    std::vector<G4double> LXe_PP    = { 7.0*eV , 7.07*eV, 7.14*eV };
    std::vector<G4double> LXe_SCINT = { 0.1, 1.0, 0.1 };
    std::vector<G4double> LXe_RIND  = { 1.59 , 1.57, 1.54 };
    std::vector<G4double> LXe_ABSL  = { 35.*cm, 35.*cm, 35.*cm}; //atten length
    std::vector<G4double> LXe_scint_e = { 0.*MeV, 10.*MeV };
    std::vector<G4double> LXe_scint_default = { 0., 120000.};
    std::vector<G4double> LXe_scint_alpha = { 0., 132000.};
    std::vector<G4double> LXe_scint_ion = { 0., 24000.};
    G4MaterialPropertiesTable *LXe_mt = new G4MaterialPropertiesTable();
    LXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", LXe_PP, LXe_SCINT);
    LXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", LXe_PP, LXe_SCINT);
    LXe_mt->AddProperty("RINDEX",        LXe_PP, LXe_RIND);
    LXe_mt->AddProperty("ABSLENGTH",     LXe_PP, LXe_ABSL);
    LXe_mt->AddProperty("ELECTRONSCINTILLATIONYIELD", LXe_scint_e, LXe_scint_default); // include QE 20%
    // and 13eV creation energy for photons - may be 15eV?
    // Fano factor assumed 1; should be much less for Xe ~ 0.13
    // but the Fano factor is already partially included in the correlated
    // electron production - therefore not the absolute Fano factor here:
    LXe_mt->AddConstProperty("ELECTRONSCINTILLATIONYIELD1", 0.);
    LXe_mt->AddConstProperty("ELECTRONSCINTILLATIONYIELD2", 1.);

    LXe_mt->AddProperty("ALPHASCINTILLATIONYIELD", LXe_scint_e, LXe_scint_alpha); // include QE 20%
    LXe_mt->AddConstProperty("ALPHASCINTILLATIONYIELD1", 1.);
    LXe_mt->AddConstProperty("ALPHASCINTILLATIONYIELD2", 0.);

    LXe_mt->AddProperty("IONSCINTILLATIONYIELD", LXe_scint_e, LXe_scint_ion); // include QE 20%
    LXe_mt->AddConstProperty("IONSCINTILLATIONYIELD1", 1.);
    LXe_mt->AddConstProperty("IONSCINTILLATIONYIELD2", 0.);

    LXe_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
    LXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",20.*ns);
    LXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",45.*ns);
    LXe->SetMaterialPropertiesTable(LXe_mt);

    std::vector<G4double> GXe_PP = { 7.0*eV, 7.07*eV, 7.14*eV };
    std::vector<G4double> GXe_SCINT = { 0.1, 1.0, 0.1 };
    std::vector<G4double> GXe_RIND  = { 1.00, 1.00, 1.00 };
    std::vector<G4double> GXe_ABSL  = { 100*m, 100*m, 100*m}; //atten length
    std::vector<G4double> GXe_scint_e = { 0.*MeV, 10.*MeV };
    std::vector<G4double> GXe_scint_default = { 0., 120000.};
    std::vector<G4double> GXe_scint_alpha = { 0., 132000.};
    std::vector<G4double> GXe_scint_ion = { 0., 24000.};
    G4MaterialPropertiesTable *GXe_mt = new G4MaterialPropertiesTable();
    GXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", GXe_PP, GXe_SCINT);
    GXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", GXe_PP, GXe_SCINT);
    GXe_mt->AddProperty("RINDEX",        GXe_PP, GXe_RIND);
    GXe_mt->AddProperty("ABSLENGTH",     GXe_PP, GXe_ABSL);
    GXe_mt->AddConstProperty("SCINTILLATIONYIELD",12000./MeV); // include QE 20%

    GXe_mt->AddProperty("ELECTRONSCINTILLATIONYIELD", GXe_scint_e, GXe_scint_default); // include QE 20%
    GXe_mt->AddConstProperty("ELECTRONSCINTILLATIONYIELD1", 0.);
    GXe_mt->AddConstProperty("ELECTRONSCINTILLATIONYIELD2", 1.);

    GXe_mt->AddProperty("ALPHASCINTILLATIONYIELD", GXe_scint_e, GXe_scint_alpha); // include QE 20%
    GXe_mt->AddConstProperty("ALPHASCINTILLATIONYIELD1", 1.);
    GXe_mt->AddConstProperty("ALPHASCINTILLATIONYIELD2", 0.);

    GXe_mt->AddProperty("IONSCINTILLATIONYIELD", LXe_scint_e, LXe_scint_ion); // include QE 20%
    GXe_mt->AddConstProperty("IONSCINTILLATIONYIELD1", 1.);
    GXe_mt->AddConstProperty("IONSCINTILLATIONYIELD2", 0.);

    GXe_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
    GXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",20.*ns);
    GXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",45.*ns);
    GXe->SetMaterialPropertiesTable(GXe_mt);


    lXeMat = nist->FindOrBuildMaterial("LXe");
    gXeMat = nist->FindOrBuildMaterial("GXe");

    // making quartz
    G4Element* O  = new G4Element
        (name="Oxygen"  ,symbol="O" , z= 8., a=16.00*g/mole);
    G4Element* Si = new G4Element
        (name="Silicon",symbol="Si" , z= 14., a=28.09*g/mole);
    G4Material* quartz = new G4Material
        (name="quartz", density=2.200*g/cm3, ncomponents=2);
    quartz->AddElement(Si, 1);
    quartz->AddElement(O , 2);

    G4MaterialPropertiesTable *quartz_mt = new G4MaterialPropertiesTable();
    quartz_mt->AddProperty("RINDEX", { 5.0*eV, 6.69*eV, 7.50*eV }, { 1.51, 1.57, 1.61 });
    quartz_mt->AddProperty("ABSLENGTH", { 5.0*eV, 6.69*eV, 7.50*eV }, { 3.0*cm, 3.0*cm, 3.0*cm });
    quartz->SetMaterialPropertiesTable(quartz_mt);

    quartzMat = nist->FindOrBuildMaterial("quartz");
    

    G4Material* vacuum = new G4Material 
    (name="Vacuum", z=1., a=1.*g/mole, density=1.e-20*g/cm3,
     kStateGas, temperature=0.1*kelvin, pressure=1.e-20*bar);
     worldMat=nist->FindOrBuildMaterial("Vacuum");

        // aluminium
    G4Element* Al = new G4Element
        (name="Aluminium"  ,symbol="Al" , z= 13., a=26.98*g/mole);  
    G4Material* metalAl = new G4Material
        (name="MetalAluminium", density=2.700*g/cm3, ncomponents=1);
    metalAl->AddElement(Al, 1);


    // photocathode aluminium
    G4Material* cathmetalAl = new G4Material
        (name="CathodeMetalAluminium", density=2.700*g/cm3, ncomponents=1);
    cathmetalAl->AddElement(Al, 1);


    G4MaterialPropertiesTable* qeTable = LoadQuantumEfficiency("eff.dat");
    if (qeTable) {
        cathmetalAl->SetMaterialPropertiesTable(qeTable);
    }
    photocathodeMat = nist->FindOrBuildMaterial("CathodeMetalAluminium");
    
    // PTFE
    // Define elements
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elF = nist->FindOrBuildElement("F");
    // Define PTFE material (C2F4)n
    G4Material* PTFE = new G4Material("PTFE", 2.2*g/cm3, 2);  // Density ~2.2 g/cm3
    PTFE->AddElement(elC, 2);
    PTFE->AddElement(elF, 4);
    G4MaterialPropertiesTable* ptfeMPT = new G4MaterialPropertiesTable();
    const G4int nEntries = 2;
    G4double photonEnergy[nEntries] = {1.0*eV, 7.5*eV};  // Covers all energies here
    G4double reflectivity[nEntries] = {1.0, 1.0};  // 100% reflective
    ptfeMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, nEntries);

    PTFE->SetMaterialPropertiesTable(ptfeMPT);

    PTFEMat = nist->FindOrBuildMaterial("PTFE");

      // titanium
    G4double z = 22.;
    G4double a = 47.867 * g/mole;
    G4double density = 4.54 * g/cm3;
    
    G4Material* titanium = new G4Material("Titanium", z, a, density);

    cryoMat = nist->FindOrBuildMaterial("Titanium");

    // copper
    G4Element* Cu = new G4Element
        (name="Copper"  ,symbol="Cu" , z= 29., a=63.55*g/mole);  
    G4Material* metalCu = new G4Material
        (name="MetalCopper", density=8.960*g/cm3, ncomponents=1);
    metalCu->AddElement(Cu, 1);
    copperMat = nist->FindOrBuildMaterial("MetalCopper");

    // kovar
    // Define elements
    G4Element* elFe = nist->FindOrBuildElement("Fe");
    G4Element* elNi = nist->FindOrBuildElement("Ni");
    G4Element* elCo = nist->FindOrBuildElement("Co");
    G4Element* elMn = nist->FindOrBuildElement("Mn");

    G4double kovarDensity = 8.36 * g/cm3;
    G4Material* kovar = new G4Material("Kovar", kovarDensity, 4);
    kovar->AddElement(elFe, 53.7 * perCent);
    kovar->AddElement(elNi, 29.0 * perCent);
    kovar->AddElement(elCo, 17.0 * perCent);
    kovar->AddElement(elMn, 0.3 * perCent);
    kovarMat = nist->FindOrBuildMaterial("Kovar");


}

MyDetectorConstruction::~MyDetectorConstruction()
{
}

void MyDetectorConstruction::ConstructSDandField() {

    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");

    sdManager->AddNewDetector(sensDet);  // Register it
    if (logicPhotocathode) {
        logicPhotocathode->SetSensitiveDetector(sensDet);
        G4cout << "Sensitive detector attached to Photocathode volume!" << G4endl;
    } else {
        G4cerr << "Error: logicPhotocathode not found!" << G4endl;
    }
}

    G4VPhysicalVolume* MyDetectorConstruction::Construct()
    {
    // Initialize vis attributes
    G4Colour glassyBlue(0.0, 0.9, 0.9, 1); // R, G, B, Alpha
    G4VisAttributes* glassVis = new G4VisAttributes(glassyBlue);
    glassVis->SetForceSolid(true);
    G4Colour white(1., 1., 1., 1); // R, G, B, Alpha
    G4VisAttributes* whiteVis = new G4VisAttributes(white);
    whiteVis->SetForceSolid(true);
    G4Colour red(1., 0., 0., 0.8); // R, G, B, Alpha
    G4VisAttributes* redVis = new G4VisAttributes(red);
    redVis->SetForceSolid(true);

    //Initialize world volume
    G4double worldSize = 1.0*m;
    solidWorld = new G4Box("solidWorld", worldSize, worldSize, worldSize);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "physWorld", 0, false, 0, true);

    // TPC Params
    G4int numSides = 12; // Dodecagon
    G4double xenonRadius = 0.25*m;
    G4double gxeHeight = 0.055*m;
    G4double lxeHeight = 0.524*m;
    G4double xenonHeight = lxeHeight + gxeHeight ; //+ 3.7 * cm;
    // NEED TO MAKE THE BIT BETWEEN PMTS AND XENON LET LIGHT TRAVEL


    // Classic medium purple, semi-opaque
    // G4Colour purpleColor(0.6, 0.0, 0.8, 0.8); // R, G, B, Alpha
    // G4VisAttributes* purpleVis = new G4VisAttributes(purpleColor);
    // purpleVis->SetForceSolid(true);
    // G4Colour lilacColor(0.5, 0.0, 1., 0.8); // R, G, B, Alpha
    // G4VisAttributes* lilacVis = new G4VisAttributes(lilacColor);
    // lilacVis->SetForceSolid(true);


    // Wall Params
    G4double wallThickness = 0.0285*m;
    G4double wallOuterRadius = xenonRadius + wallThickness;

    // PMT parameters
    G4double pmtRadius = 0.028 * m;
    G4double pmtThickness = 0.005 * m;
    G4double photocathodeRadius = 22.5 * mm;
    G4double photocathodeThickness = 0.002 * m;
    G4double kovarHeight = pmtThickness + photocathodeThickness + 0.001*cm;

    // PTFE Tile On Top/Bottom Params
    // PTFE tile shape: thin dodecagon plate
    G4double ptfeTileThickness = 0.005 * m;
    G4double tileInnerRadius = 0.0;             // center filled
    G4double tileOuterRadius = wallOuterRadius; // slightly beyond active volume
    G4double startAngle = 0.*deg;
    G4double totalAngle = 360.*deg;

    // Field rings
    G4double ringThickness = 12.7 * mm;             // Radial thickness of ring
    G4double ringRadius = xenonRadius + wallThickness / 2;  // in the wall
    G4double ringHeight = 3.2 * mm;                // Axial height of each ring
    G4double activeHeight = xenonHeight;
    G4int numRings = 48;  // Number of rings
    // Spacing between rings - evenly spaced
    G4double ringSpacing = 1 * cm;


    // Inner Cryo Params
    G4double innerCryoOuterRadius = (24.25 * 2.54 * cm) / 2;   // 24.25 inch diameter
    G4double innerCryoHeight = 39.75 * 2.54 * cm;              // 39.75 inch tall
    G4ThreeVector cavityPos = G4ThreeVector(0., 0., 0.);

    // Vaccuum between cryostats params
    G4double vacuumMargin = 20.0 * cm;  // Guess: can adjust buffer space

    G4double vacuumRadius = innerCryoOuterRadius + vacuumMargin;
    G4double vacuumHeight = innerCryoHeight + 2 * 20 * cm;

    // Outer Cryo Params
    G4double cryostatThickness = 0.0056642 * m;  // Cryostat thickness
    G4double outerCryoRadius = vacuumRadius + cryostatThickness; //Guess
    G4double outerCryoHeight = innerCryoHeight + 2 * 20 * cm + 2 * cryostatThickness;

    // Copper params
    G4double copperDiskRadius = 27.5 * cm;  // 55 cm diameter → 27.5 cm radius
    G4double copperDiskThickness = 5.0 * cm;


    // ----------------------------------------------------------
    // --- TPC ---
    // Define dodecagon xenon volume (12-sided prism)
    G4double zPlanes[2] = { -xenonHeight/2, xenonHeight/2 };
    G4double rInner[2] = { 0.0, 0.0 };
    G4double rOuter[2] = { xenonRadius, xenonRadius };

    G4Polyhedra* solidXenon = new G4Polyhedra("solidXenon", 0.*deg, 360.*deg, numSides, 2, zPlanes, rInner, rOuter);
    G4LogicalVolume* logicXenon = new G4LogicalVolume(solidXenon, worldMat, "logicXenon");
    G4VPhysicalVolume* physXenon = new G4PVPlacement(nullptr, {}, logicXenon, "physXenon", logicWorld, false, 0, true);

    // Liquid Xenon (lower part)
    G4double lxeZPlanes[2] = { -lxeHeight/2, lxeHeight/2 };
    G4double lxeROuter[2] = { xenonRadius, xenonRadius };
    G4Polyhedra* solidLXe = new G4Polyhedra("solidLXe", 0.*deg, 360.*deg, numSides, 2, lxeZPlanes, rInner, lxeROuter);
    G4LogicalVolume* logicLXe = new G4LogicalVolume(solidLXe, lXeMat, "logicLXe");
    //logicLXe->SetVisAttributes(purpleVis);

    G4ThreeVector lxePos(0., 0., -gxeHeight/2);
    new G4PVPlacement(nullptr, lxePos, logicLXe, "physLXe", logicXenon, false, 0, true);

    // Gaseous Xenon (upper part)
    G4double gxeZPlanes[2] = { -gxeHeight/2, gxeHeight/2 };
    G4double gxeROuter[2] = { xenonRadius, xenonRadius };
    G4Polyhedra* solidGXe = new G4Polyhedra("solidGXe", 0.*deg, 360.*deg, numSides, 2, gxeZPlanes, rInner, gxeROuter);
    G4LogicalVolume* logicGXe = new G4LogicalVolume(solidGXe, gXeMat, "logicGXe");
    //logicGXe->SetVisAttributes(lilacVis);

    G4ThreeVector gxePos(0., 0., lxeHeight/2);
    new G4PVPlacement(nullptr, gxePos, logicGXe, "physGXe", logicXenon, false, 0, true);

    fScoringVolume = logicLXe;

    // --- Field Rings ---
    // Define the ring solid
    G4Tubs* solidRing = new G4Tubs("solidRing", ringRadius, ringRadius + ringThickness, ringHeight / 2, 0., 360.*deg);

    // Create logical volume with copper material
    G4LogicalVolume* logicRing = new G4LogicalVolume(solidRing, copperMat, "logicRing");

    // Set copper colour
    // G4Colour copperColor(0.72, 0.45, 0.20, 0.6); // R, G, B, Alpha
    // G4VisAttributes* copperVis = new G4VisAttributes(copperColor);
    // copperVis->SetForceSolid(false);
    // logicRing->SetVisAttributes(copperVis);

    // --- Place the field rings ---
    for (G4int i = 1; i < numRings+1; ++i) {
        // Centered on Z = 0; distribute along the Z-axis
        G4double zPos = (-activeHeight / 2) + i * ringSpacing;
        //G4cout << "zPos of ring " << i - 1 << " is: " << zPos << G4endl;


        // Place each field ring in the world
        new G4PVPlacement(
            nullptr,
            G4ThreeVector(0., 0., zPos),        // Position of the ring along Z-axis
            logicRing,                           // Logical volume for the ring
            "physRing",                          // Name of the physical volume
            logicWorld,                          // The world logical volume
            false,
            i,                                   // Unique ID for each ring
            true                                 // Check overlap
        );

    }

    // --- PTFE wall --- (around xenon, with holes for rings)
    G4double wallROuter[2] = { wallOuterRadius, wallOuterRadius };
    G4double wallZPlanes[2] = { -xenonHeight/2 , xenonHeight/2};

    // Polyhedron for the PTFE wall (with no holes yet)
    G4Polyhedra* solidWall = new G4Polyhedra("solidWall", 0.*deg, 360.*deg, numSides, 2, wallZPlanes, rOuter, wallROuter);

    // Create a copy of the wall to modify (subtract rings)
    G4VSolid* wallWithHoles = solidWall;

    // Loop over all rings and subtract them from the wall
    for (G4int i = 1; i < numRings+1; ++i) {
        // Centered on Z = 0; distribute along the Z-axis
        G4double zPos = (-activeHeight / 2) + i * ringSpacing ;

        // Position of the current ring
        G4ThreeVector ringPos(0., 0., zPos);

        // Subtract the current ring from the wall
        wallWithHoles = new G4SubtractionSolid("wallWithHoles", wallWithHoles, solidRing, nullptr, ringPos);
    }

    // Logical volume for the PTFE wall (with holes)
    G4LogicalVolume* logicWall = new G4LogicalVolume(wallWithHoles, PTFEMat, "logicWall");

    // Place the PTFE wall with holes (around the Xenon volume)
    new G4PVPlacement(nullptr, {}, logicWall, "physWall", logicWorld, false, 0, true);


    // --- PMTs ---
    solidPMT = new G4Tubs("solidPMT", 0., pmtRadius, pmtThickness/2, 0., 360.*deg);
    logicPMT = new G4LogicalVolume(solidPMT, quartzMat, "logicPMT");
    //logicPMT->SetVisAttributes(redVis);

    // Photocathodes
    // Define the solid photocathode
    G4Tubs* solidPhotocathode = new G4Tubs("solidPhotocathode", 0., photocathodeRadius, photocathodeThickness/2, 0., 360.*deg);
    logicPhotocathode = new G4LogicalVolume(solidPhotocathode, photocathodeMat, "logicPhotocathode");
    //logicPhotocathode->SetVisAttributes(whiteVis);

    // Define the solid kovar (outer casing)
    G4Tubs* solidKovar = new G4Tubs("solidKovar", 0., pmtRadius, kovarHeight/2, 0., 360.*deg);

    // Position of the subtraction solid (cut the cavity from the top)
    G4ThreeVector cavityPosCasing(0, 0, (kovarHeight - photocathodeThickness) / 2);

    // Subtract the cavity from the top of the kovar to make room for the photocathode
    G4SubtractionSolid* solidCasing = new G4SubtractionSolid("solidCasing", solidKovar, solidPhotocathode, nullptr, cavityPosCasing);
    G4LogicalVolume* logicCasing = new G4LogicalVolume(solidCasing, kovarMat, "logicCasing");
    //logicCasing->SetVisAttributes(glassVis);

    // Define non-overlapping spacing
    G4double pitch = 2.4 * pmtRadius;              // No overlap
    G4double dx = pitch * std::sqrt(3)/2.0;        // Horizontal spacing in hex grid
    G4double dy = pitch * 0.75;                    // Vertical spacing

    int gridRadius = 4;  // 4 gives 5 rings = 61 PMTs (not sure why -1)
    int pmtIndex = 0;

    // Define z positions
    G4double assemblySpacing = xenonHeight + 2*pmtThickness; // Distance between the two assemblies

    G4double casingZ = (kovarHeight + assemblySpacing) / 2.0;
    G4double cathodeZ = casingZ - (kovarHeight - photocathodeThickness)/2;  // match subtraction


    for (int q = -gridRadius; q <= gridRadius; ++q) {
         int r1 = std::max(-gridRadius, -q - gridRadius);
         int r2 = std::min(gridRadius, -q + gridRadius);

        for (int r = r1; r <= r2; ++r) {
            // Convert axial to cartesian
            G4double x = dx * (q + r/2.0);
            G4double y = dy * r;

            // Save position for later PTFE hole subtraction
            pmtPositions.emplace_back(x, y, 0.);

            // --- Top PMT ---
            // Place top assembly (flipped)
            G4RotationMatrix* rot180 = new G4RotationMatrix();
            rot180->rotateX(180.*deg);

            new G4PVPlacement(rot180,
                            G4ThreeVector(x,y,casingZ),   // mirror placement
                            logicCasing, "CasingUp",
                            logicWorld, false, pmtIndex, true);

            new G4PVPlacement(rot180,
                            G4ThreeVector(x,y,cathodeZ),  // match the flipped cavity
                            logicPhotocathode, "PhotocathodeUp",
                            logicWorld, false, pmtIndex, true);

            physPMT_Top = new G4PVPlacement(rot180,
                G4ThreeVector(x, y, xenonHeight/2 + pmtThickness/2),
                logicPMT, "physPMT_Top", logicWorld, false, pmtIndex, true);

            // --- Bottom PMT ---
            // Place bottom assembly (upright)
            new G4PVPlacement(nullptr,                      // no rotation
                G4ThreeVector(x,y,-casingZ),   // position
                logicCasing, "CasingDown",     // logic and name
                logicWorld, false, pmtIndex+1000, true);

            new G4PVPlacement(nullptr,
                G4ThreeVector(x,y,-cathodeZ),  // position matches cavity
                logicPhotocathode, "PhotocathodeDown",
                logicWorld, false, pmtIndex+1000, true);

            physPMT_Bottom = new G4PVPlacement(nullptr,
                G4ThreeVector(x, y, -xenonHeight/2 - pmtThickness/2),
                logicPMT, "physPMT_Bottom", logicWorld, false, pmtIndex + 1000, true);

            ++pmtIndex;
        }
    }

    //  --- Make PTFE around the PMTs ---
    // This creates a thin 12-sided polygonal tile
    G4Polyhedra* baseTile = new G4Polyhedra("baseTile",
                                            startAngle,
                                            totalAngle,
                                            numSides,
                                            2,  // z-section count
                                            (G4double[]){ -ptfeTileThickness/2, ptfeTileThickness/2 },  // z positions
                                            (G4double[]){ tileInnerRadius, tileInnerRadius },           // inner radii
                                            (G4double[]){ tileOuterRadius, tileOuterRadius });          // outer radii

    // Start with the base tile
    G4VSolid* tileWithHoles = baseTile;

    // Create a PMT hole solid to subtract 
    G4double holeRadius = pmtRadius;
    G4Tubs* pmtHole = new G4Tubs("pmtHole", 0., holeRadius, ptfeTileThickness, 0., 360.*deg);

    for (const auto& pos : pmtPositions) {
        G4ThreeVector holePos = G4ThreeVector(pos.x(), pos.y(), 0.);  // holes in XY plane

        tileWithHoles = new G4SubtractionSolid("tileWithHoles", tileWithHoles, pmtHole,
                                            nullptr, holePos);
    }

    // Logical volume from final solid
    G4LogicalVolume* logicPTFETile = new G4LogicalVolume(tileWithHoles, PTFEMat, "logicPTFETile");

    //logicPTFETile->SetVisAttributes(glassVis);

    // Top tile placement
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0., 0., xenonHeight/2 + ptfeTileThickness/2),
        logicPTFETile,
        "PTFE_Tile_Top",
        logicWorld,
        false,
        0,
        true
    );

    // Bottom tile placement
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0., 0., -xenonHeight/2 - ptfeTileThickness/2),
        logicPTFETile,
        "PTFE_Tile_Bottom",
        logicWorld,
        false,
        1,
        true
    );

    // --- Copper Blocks ---
    // Define copper disk solid
    G4Tubs* solidCopperDisk = new G4Tubs("solidCopperDisk", 0., copperDiskRadius, copperDiskThickness / 2, 0., 360. * deg);

    // Logical volume for copper disk
    G4LogicalVolume* logicCopperDisk = new G4LogicalVolume(solidCopperDisk, copperMat, "logicCopperDisk");

    //logicCopperDisk->SetVisAttributes(copperVis);


    // --- Top copper disk (above PMTs and PTFE tile) ---
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0., 0., xenonHeight/2 + pmtThickness + kovarHeight + copperDiskThickness / 2),
        logicCopperDisk,
        "CopperDisk_Top",
        logicWorld,
        false,
        0,
        true
    );

    // Bottom copper disk (below PMTs and PTFE tile)
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0., 0., -xenonHeight/2 - pmtThickness - kovarHeight - copperDiskThickness / 2),
        logicCopperDisk,
        "CopperDisk_Bottom",
        logicWorld,
        false,
        1,
        true
    );
    // The inner cut-out where detector goes - TODO make this line up with the PTFE
    G4double cavityRadius = innerCryoOuterRadius - 0.0056642 * m;       //
    G4double cavityHeight = xenonHeight + 2 * pmtThickness + 2.0 * kovarHeight + 2 * copperDiskThickness; // sits on top of copper disk


    // --- LXe between PTFE Wall and Inner Cryostat ---
    G4Tubs* solidLXeBlock = new G4Tubs("solidLXeBlock", 0., innerCryoOuterRadius - 0.0056642 * m, xenonHeight / 2, 0., 360.*deg);
    G4Polyhedra* cavityXenon = new G4Polyhedra("cavityXenon", 0.*deg, 360.*deg, numSides, 2, wallZPlanes, rInner, wallROuter);

    G4SubtractionSolid* solidOuterLXe = new G4SubtractionSolid("solidOuterLXe", solidLXeBlock, cavityXenon, nullptr, cavityPos);

    G4LogicalVolume* logicOuterLXe = new G4LogicalVolume(solidOuterLXe, lXeMat, "logicOuterLXe");

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicOuterLXe, "physOuterLXe", logicWorld, false, 0, true);

    // --- Inner Cryostat ---
    // Overall solid block size (inner cryostat)
    G4Tubs* solidCryoBlock1 = new G4Tubs("solidCryoBlock1", 0., innerCryoOuterRadius, innerCryoHeight / 2, 0., 360.*deg);

    // Inner cavity (cylinder to subtract)
    G4Tubs* solidCavity = new G4Tubs("solidCavity", 0., cavityRadius, cavityHeight / 2, 0., 360.*deg);

    // Subtract the cavity from the block
    // Make sure the cavity is centered inside the block (adjust if needed)
    G4SubtractionSolid* solidInnerCryo = new G4SubtractionSolid("solidInnerCryo", solidCryoBlock1, solidCavity, nullptr, cavityPos);

    // Create logical volume with steel material
    G4LogicalVolume* logicInnerCryo = new G4LogicalVolume(solidInnerCryo, cryoMat, "logicInnerCryo");
    //logicInnerCryo->SetVisAttributes(grayVis);

    // Place it into the world
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicInnerCryo, "physInnerCryo", logicWorld, false, 0, true);

    // --- Vacuum Between Cryostats ---

    // Define solid cylinder of vacuum
    G4Tubs* solidVacuumBlock = new G4Tubs("solidVacuumBlock", 0., vacuumRadius, vacuumHeight / 2, 0., 360.*deg);

    G4SubtractionSolid* solidVacuum = new G4SubtractionSolid("solidVacuum", solidVacuumBlock, solidCryoBlock1, nullptr, cavityPos);

    G4LogicalVolume* logicVacuum = new G4LogicalVolume(solidVacuum, worldMat, "logicVacuum");

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicVacuum, "physVacuum", logicWorld, false, 0, true);

    // --- Outer Cryostat ---
    // Define solid cylinder of vacuum
    G4Tubs* solidCryoBlock2 = new G4Tubs("solidCryoBlock2", 0., outerCryoRadius, outerCryoHeight / 2, 0., 360.*deg);

    G4SubtractionSolid* solidOuterCryo = new G4SubtractionSolid("solidOuterCryo", solidCryoBlock2, solidVacuumBlock, nullptr, cavityPos);

    G4LogicalVolume* logicOuterCryo= new G4LogicalVolume(solidOuterCryo, cryoMat, "logicOuterCryo");
    //logicOuterCryo->SetVisAttributes(grayVis);

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicOuterCryo, "physOuterCryo", logicWorld, false, 0, true);


    return physWorld;
}

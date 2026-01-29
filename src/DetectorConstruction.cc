#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"

#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
//#include "CADMesh.hh"
#include "G4Element.hh"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

DetectorConstruction* DetectorConstruction::fInstance = nullptr;

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
    fMessenger(nullptr),
    fCustomRockDensity(2.7 * g / cm3),
    fUseCustomRock(false),
    fLabName("default"),
    fFilePattern("")
{
    fMessenger = new DetectorMessenger(this);
    fInstance = this;
}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
    if (fInstance == this) {
        fInstance = nullptr;
    }
}

void DetectorConstruction::SetRockDensity(G4double density)
{
    fCustomRockDensity = density;
    G4cout << "==================================================" << G4endl;
    G4cout << "Rock density set to: " << fCustomRockDensity / (g / cm3) << " g/cm3" << G4endl;
    G4cout << "==================================================" << G4endl;
}

void DetectorConstruction::SetRockComposition(const G4String& filename)
{   
    //G4String rockCompDir = "./Rock Composition Data/";
    G4String rockCompDir = "../../Rock Composition Data/";
    LoadRockComposition(rockCompDir + filename);
    fUseCustomRock = true;
    PrintRockComposition();
}

void DetectorConstruction::SetLabName(const G4String& name)
{
    fLabName = name;
    G4cout << "==================================================" << G4endl;
    G4cout << "Lab name set to: " << fLabName << G4endl;
    G4cout << "==================================================" << G4endl;
}

void DetectorConstruction::LoadRockComposition(const G4String& filename)
{
    fRockComposition.clear();

    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        G4cerr << "ERROR: Could not open rock composition file: " << filename << G4endl;
        return;
    }

    std::string line;
    G4int lineNum = 0;
    G4double totalFraction = 0.0;

    while (std::getline(inFile, line)) {
        lineNum++;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string element;
        G4double fraction;

        if (iss >> element >> fraction) {
            fRockComposition[element] = fraction;
            totalFraction += fraction;
        }
        else {
            G4cerr << "WARNING: Could not parse line " << lineNum << ": " << line << G4endl;
        }
    }

    inFile.close();

    G4cout << "Loaded rock composition with " << fRockComposition.size()
        << " elements from file: " << filename << G4endl;
    G4cout << "Total mass fraction: " << totalFraction << G4endl;

    if (std::abs(totalFraction - 1.0) > 0.01) {
        G4cout << "WARNING: Mass fractions do not sum to 1.0!" << G4endl;
    }
}

void DetectorConstruction::PrintRockComposition()
{
    G4cout << "==================================================" << G4endl;
    G4cout << "Custom Rock Composition:" << G4endl;
    G4cout << "==================================================" << G4endl;
    G4cout << "Density: " << fCustomRockDensity / (g / cm3) << " g/cm3" << G4endl;
    G4cout << "--------------------------------------------------" << G4endl;
    G4cout << std::setw(10) << "Element" << std::setw(20) << "Mass Fraction" << G4endl;
    G4cout << "--------------------------------------------------" << G4endl;

    for (const auto& pair : fRockComposition) {
        G4cout << std::setw(10) << pair.first
            << std::setw(20) << std::fixed << std::setprecision(6)
            << pair.second << G4endl;
    }
    G4cout << "==================================================" << G4endl;
}

G4Material* DetectorConstruction::CreateCustomRock()
{
    if (fRockComposition.empty()) {
        G4cerr << "ERROR: No rock composition loaded!" << G4endl;
        return nullptr;
    }

    G4NistManager* nist = G4NistManager::Instance();

    G4Material* customRock = new G4Material("CustomRock", fCustomRockDensity,
        fRockComposition.size(), kStateSolid);

    for (const auto& pair : fRockComposition) {
        G4Element* element = nist->FindOrBuildElement(pair.first);
        if (element) {
            customRock->AddElement(element, pair.second);
        }
        else {
            G4cerr << "ERROR: Could not find element: " << pair.first << G4endl;
        }
    }

    G4cout << "Created custom rock material with " << fRockComposition.size()
        << " elements" << G4endl;

    return customRock;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();


    // Defining the rock composition for the Hideout and Quantum Lab

    G4Element* Si = nist->FindOrBuildElement("Si");
    G4Element* Ca = nist->FindOrBuildElement("Ca");
    G4Element* O = nist->FindOrBuildElement("O");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* Al = nist->FindOrBuildElement("Al");
    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* K = nist->FindOrBuildElement("K");
    G4Element* Na = nist->FindOrBuildElement("Na");
    G4Element* Mg = nist->FindOrBuildElement("Mg");
    G4Element* Fe = nist->FindOrBuildElement("Fe");
    G4Element* P = nist->FindOrBuildElement("P");
    G4Element* S = nist->FindOrBuildElement("S");
    G4Element* F = nist->FindOrBuildElement("F");
    G4Element* Cl = nist->FindOrBuildElement("Cl");


    G4double hideoutDensity = 2.6949673 * g / cm3;  // hideout rock density
    G4Material* hideoutRock = new G4Material("hideoutRock", hideoutDensity, 12, kStateSolid);  // 12 components

    // Add elements with correct fractions
    hideoutRock->AddElement(Si, 0.362176);  // silicon
    hideoutRock->AddElement(Ca, 0.002794);  // calcium
    hideoutRock->AddElement(O, 0.498402891);  // oxygen
    hideoutRock->AddElement(C, 0.000417109);  // carbon
    hideoutRock->AddElement(Al, 0.058924);  // aluminum
    hideoutRock->AddElement(H, 0.000958);  // hydrogen
    hideoutRock->AddElement(K, 0.0165);  // potassium
    hideoutRock->AddElement(Na, 0.034698);  // sodium
    hideoutRock->AddElement(Mg, 0.010806);  // magnesium
    hideoutRock->AddElement(Fe, 0.01237);  // iron
    hideoutRock->AddElement(S, 0.001143);  // sulfur
    hideoutRock->AddElement(F, 0.000811);  // fluorine

    // vacuum for the space between particle start position and rock
    G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");


    // inner and outer radius of sphere for world volume
    G4double innerRadius = 0.0 * m;
    G4double outerRadius = 8.4 * m;

    // Define start and end angles for theta and phi to create a hemisphere...note deltaTheta / 2 gives the upper hemisphere
    G4double startPhi = 0.0 * deg;
    G4double deltaPhi = 360.0 * deg;
    G4double startTheta = 0.0 * deg;
    G4double deltaTheta = 180.0 * deg;
    // Create the world spherical shape
    G4Sphere* worldS = new G4Sphere("WorldSphere", innerRadius, outerRadius, startPhi, deltaPhi, startTheta, deltaTheta / 2);
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, vacuum, "WorldSphere"); // where we assign the vacuum to world volume
    G4VPhysicalVolume* worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "WorldSphere", 0, false, 0, true);


    // Colors
    G4VisAttributes* visRed = new G4VisAttributes(G4Colour::Red());
    G4VisAttributes* visBlue = new G4VisAttributes(G4Colour::Blue());
    G4VisAttributes* visYellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    G4VisAttributes* visTrans = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.0));
    //G4VisAttributes* visWhite = new G4VisAttributes(G4Colour::White());
    G4VisAttributes* visWhite = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.08));
    G4VisAttributes* visGreen = new G4VisAttributes(G4Colour::Green());

    G4VisAttributes* visWorld = new G4VisAttributes(G4Colour::White());
    worldLV->SetVisAttributes(visBlue);
    //visBlue->SetForceSolid(true);




    G4double labX = 3.049999 / 2. * m;
    G4double labY = 5.199999 / 2. * m;
    G4double labZ = 3.049999 / 2. * m;

    G4VSolid* b3S = new G4Box("lab", labX, labY, labZ);
    G4LogicalVolume* boxLV = new G4LogicalVolume(b3S, air, "labLV");
    G4VPhysicalVolume* boxPV = new G4PVPlacement(0, G4ThreeVector(0, 0, labZ), boxLV, "labPV", worldLV, false, 0, true);
    boxLV->SetVisAttributes(visBlue);

    // Dimensions for the rectangular shell
    G4double innerLengthX = 3.05 * m;
    G4double innerLengthY = 5.2 * m;
    G4double innerLengthZ = 3.05 * m;

    G4double thicknessTop = 5.0 * m;     // Thickness above the lab space
    G4double thicknessSide = 5.0 * m;    // Thickness on sides
    G4double thicknessBelow = 0.0 * m;   // No thickness below the lab space

    // Outer box dimensions
    G4double outerLengthX = innerLengthX + 2.0 * thicknessSide;
    G4double outerLengthY = innerLengthY + 2.0 * thicknessSide;
    G4double outerLengthZ = innerLengthZ + thicknessTop; // Only add thickness to top

    G4double innerRadius2 = 0.0 * m;
    G4double outerRadius2 = 8.3 * m;

    // ========== CUSTOM ROCK IMPLEMENTATION ==========
    // Determine which rock material to use
    G4Material* rockMaterial;
    if (fUseCustomRock) {
        rockMaterial = CreateCustomRock();
        if (!rockMaterial) {
            G4cout << "WARNING: Failed to create custom rock. Falling back to snolabRock." << G4endl;
            rockMaterial = hideoutRock;
        }
        else {
            G4cout << "Using custom rock material in geometry." << G4endl;
        }
    }
    else {
        rockMaterial = hideoutRock;  // Default material
    }
    // ================================================

    G4Sphere* rockSphere = new G4Sphere("rockSphere", innerRadius2, outerRadius2, startPhi, deltaPhi, startTheta, deltaTheta / 2);
    G4LogicalVolume* rockSphereLV = new G4LogicalVolume(rockSphere, rockMaterial, "rockSphere");
    G4Box* innerBox = new G4Box("InnerBox", innerLengthX / 2, innerLengthY / 2,
        innerLengthZ / 2);
    G4SubtractionSolid* rectangularShell =
        new G4SubtractionSolid("rShell", rockSphere, innerBox,
            0, G4ThreeVector(0, 0, 3.05 / 2 * m));


    // Define logical volume for the shell

    G4LogicalVolume* shellLV = new G4LogicalVolume(rectangularShell, rockMaterial, "rockLV");
    //G4RotationMatrix* rotation = new G4RotationMatrix();
    //rotation->rotateZ(55.0 * deg);

    // place the lab space on the equator
    G4VPhysicalVolume* shellPV =
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0 * m), shellLV, "rockPV", worldLV, false, 0, true);
    //G4VisAttributes* visShell = new G4VisAttributes(G4Colour::Red());
    //shellLV->SetVisAttributes(visShell);
    //visShell->SetForceSolid(true);

    G4VisAttributes* visShell = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.04)); // Red with 30% opacity
    visShell->SetForceSolid(true);  // Ensure the hemisphere appears solid
    shellLV->SetVisAttributes(visShell);



    return worldPV;
}

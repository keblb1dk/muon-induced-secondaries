#include "DetectorConstruction.hh"
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
#include "G4Element.hh"

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}

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
    G4Material* hideoutRock = new G4Material("hideoutRock", hideoutDensity, 10, kStateSolid);  // 12 components

    // Add elements with correct fractions
    hideoutRock->AddElement(Si, 0.387775);  // silicon....0.362176
    hideoutRock->AddElement(Ca, 0.002794);  // calcium
    hideoutRock->AddElement(O, 0.524001891);  // oxygen....0.498402891
    hideoutRock->AddElement(C, 0.000417109);  // carbon
    hideoutRock->AddElement(Al, 0.058924);  // aluminum
    hideoutRock->AddElement(H, 0.000958);  // hydrogen
    //hideoutRock->AddElement(K, 0.0165);  // potassium
    //hideoutRock->AddElement(Na, 0.034698);  // sodium
    hideoutRock->AddElement(Mg, 0.010806);  // magnesium
    hideoutRock->AddElement(Fe, 0.01237);  // iron
    hideoutRock->AddElement(S, 0.001143);  // sulfur
    hideoutRock->AddElement(F, 0.000811);  // fluorine


    
    // vacuum for the space between particle start position and rock
    G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

    // air for lab cavity...can use vacuum as well
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    

    // inner and outer radius of sphere for world volume
    G4double innerRadius = 0.0 * m;
    G4double outerRadius = 8.4 * m;

    // Define start and end angles for theta and phi to create a hemisphere...note deltaTheta / 2 gives the upper hemisphere
    G4double startPhi = 0.0 * deg;
    G4double deltaPhi = 360.0 * deg;
    G4double startTheta = 0.0 * deg;
    G4double deltaTheta = 180.0 * deg;
    // Create the world hemispherical shape
    G4Sphere* worldS = new G4Sphere("WorldSphere", innerRadius, outerRadius, startPhi, deltaPhi, startTheta, deltaTheta/2);
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, vacuum, "WorldSphere"); // where we assign the vacuum to world volume
    G4VPhysicalVolume* worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "WorldSphere", 0, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());


    // Colors
    G4VisAttributes* visRed = new G4VisAttributes(G4Colour::Red());
    G4VisAttributes* visBlue = new G4VisAttributes(G4Colour::Blue());
    G4VisAttributes* visYellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));  // R, G, B, transparency 
    G4VisAttributes* visTrans = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.0)); 
    G4VisAttributes* visWhite = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.08));
    G4VisAttributes* visGreen = new G4VisAttributes(G4Colour::Green());


    // radii for the rock hemisphere
    G4double innerRadiusRock = 0.0 * m;
    G4double outerRadiusRock = 8.3 * m;

    // Create the rock hemispherical shape
    G4Sphere* rockHS = new G4Sphere("rockHS", innerRadiusRock, outerRadiusRock, startPhi, deltaPhi, startTheta, deltaTheta / 2);
    G4LogicalVolume* rockLV = new G4LogicalVolume(rockHS, hideoutRock, "rockLV"); 
    G4VPhysicalVolume* rockPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), rockLV, "rockPV", worldLV, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());


    G4VisAttributes* visRock = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2)); // red, 20% opaque
    visRock->SetForceSolid(true);    // show as semi-transparent solid
    visRock->SetForceAuxEdgeVisible(true); // draw edges so boundary is clear
    rockLV->SetVisAttributes(visRock);
    //visShell->SetForceSolid(true);

    // create lab dimensions
    G4double labX = 3.0 / 2. * m;
    G4double labY = 5.0 / 2. * m;
    G4double labZ = 3.0 / 2. * m;

    G4VSolid* labCavity = new G4Box("lab", labX, labY, labZ);
    G4LogicalVolume* boxLV = new G4LogicalVolume(labCavity, air, "labLV");
    // place lab in rock volume
    G4VPhysicalVolume* boxPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 150. * cm), boxLV, "labPV", rockLV, false, 0, true); 
    boxLV->SetVisAttributes(visGreen);

    

    return worldPV;
}

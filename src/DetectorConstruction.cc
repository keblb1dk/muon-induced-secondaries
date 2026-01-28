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


    // Defining the rock composition for the Subatomic Particle Hideout

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
    G4Material* hideoutRock = new G4Material("hideoutRock", hideoutDensity, 10, kStateSolid);  // 10 components

    // Add elements with correct fractions
    hideoutRock->AddElement(Si, 0.387775);  // silicon
    hideoutRock->AddElement(Ca, 0.002794);  // calcium
    hideoutRock->AddElement(O, 0.524001891);  // oxygen
    hideoutRock->AddElement(C, 0.000417109);  // carbon
    hideoutRock->AddElement(Al, 0.058924);  // aluminum
    hideoutRock->AddElement(H, 0.000958);  // hydrogen
    hideoutRock->AddElement(Mg, 0.010806);  // magnesium
    hideoutRock->AddElement(Fe, 0.01237);  // iron
    hideoutRock->AddElement(S, 0.001143);  // sulfur
    hideoutRock->AddElement(F, 0.000811);  // fluorine

    // Soudan
    G4double soudanDensity = 2.85 * g / cm3;  // soudan rock density
    G4Material* soudanRock = new G4Material("soudanRock", soudanDensity, 10, kStateSolid);  // 10 components

    // Add elements with correct fractions
    soudanRock->AddElement(H, 0.00304);
    soudanRock->AddElement(C, 0.00081); 
    soudanRock->AddElement(O, 0.45467); 
    soudanRock->AddElement(Na, 0.01874);
    soudanRock->AddElement(Mg, 0.0397);
    soudanRock->AddElement(Al, 0.08042); 
    soudanRock->AddElement(Si, 0.23954); 
    soudanRock->AddElement(K, 0.00334); 
    soudanRock->AddElement(Ca, 0.06513); 
    soudanRock->AddElement(Fe, 0.09461);  
	
	// Gran Sasso
	G4double lngsDensity = 2.72 * g / cm3;  // lngs rock density
    G4Material* lngsRock = new G4Material("lngsRock", lngsDensity, 8, kStateSolid);  // 8 components
    lngsRock->AddElement(H, 0.0003);
    lngsRock->AddElement(C, 0.1217); 
    lngsRock->AddElement(O, 0.5079); 
    lngsRock->AddElement(Mg, 0.0832);
    lngsRock->AddElement(Al, 0.0063); 
    lngsRock->AddElement(Si, 0.0105);  
    lngsRock->AddElement(K, 0.001);  
    lngsRock->AddElement(Ca, 0.2691); 


    
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


    G4double labX = 2.999999 / 2. * m;
    G4double labY = 4.999999 / 2. * m;
    G4double labZ = 2.999999 / 2. * m;
	
	G4VSolid* b3S = new G4Box("lab", labX, labY, labZ);
	G4LogicalVolume* boxLV = new G4LogicalVolume(b3S, air, "labLV");
	G4VPhysicalVolume* boxPV = new G4PVPlacement(0, G4ThreeVector(0, 0, labZ), boxLV, "labPV", worldLV, false, 0, true);
	boxLV->SetVisAttributes(visBlue);
	
	G4double detX = 1.0 / 2. * m;
    G4double detY = 1.0 / 2. * m;
    G4double detZ = 1.0 / 2. * nm;
	
	G4VSolid* detector = new G4Box("detector", detX, detY, detZ);
	G4LogicalVolume* detLV = new G4LogicalVolume(detector, air, "detLV");
	G4VPhysicalVolume* detPV = new G4PVPlacement(0, G4ThreeVector(0, 0, -labZ / 2 + detZ), detLV, "detPV", boxLV, false, 0, true);
	detLV->SetVisAttributes(visRed);
	
	// Dimensions for the rectangular shell
	G4double innerLengthX = 3.0 * m;
	G4double innerLengthY = 5.0 * m;
	G4double innerLengthZ = 3.0 * m;
	
	G4double thicknessTop = 5.0 * m;     // Thickness above the lab space
	G4double thicknessSide = 5.0 * m;    // Thickness on sides
	G4double thicknessBelow = 0.0 * m;   // No thickness below the lab space
	
	// Outer box dimensions
	G4double outerLengthX = innerLengthX + 2.0 * thicknessSide;
	G4double outerLengthY = innerLengthY + 2.0 * thicknessSide;
	G4double outerLengthZ = innerLengthZ + thicknessTop; // Only add thickness to top
	G4double innerRadius2 = 0.0 * m;
	G4double outerRadius2 = 8.3 * m;
	G4Sphere* rockSphere = new G4Sphere("rockSphere", innerRadius2, outerRadius2, startPhi, deltaPhi, startTheta, deltaTheta / 2);
	G4LogicalVolume* rockSphereLV = new G4LogicalVolume(rockSphere, lngsRock, "rockSphere");
	G4Box* innerBox = new G4Box("InnerBox", innerLengthX / 2, innerLengthY / 2,
		innerLengthZ / 2);
	G4SubtractionSolid* rectangularShell =
		new G4SubtractionSolid("rShell", rockSphere, innerBox,
			0, G4ThreeVector(0, 0, 3.0 / 2 * m));
	  
	// Define logical volume for the shell
	G4LogicalVolume* shellLV = new G4LogicalVolume(rectangularShell, lngsRock, "rShell");
	
	// place the lab space on the equator
	G4VPhysicalVolume* shellPV =
		new G4PVPlacement(0, G4ThreeVector(0, 0, 0 * m), shellLV, "rShell", worldLV, false, 0, true);
		
	
	G4VisAttributes* visShell = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.04)); // Red with 30% opacity
	visShell->SetForceSolid(true);  // Ensure the hemisphere appears solid
	shellLV->SetVisAttributes(visShell);
    

    return worldPV;
}

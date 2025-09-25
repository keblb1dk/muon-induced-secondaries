#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh" 
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <iostream>

//create variables for the angular sampling
struct MuonData {
    double zenith;
    double azimuth;
    double intensity;
    double cdf;
};

//create variables for the energy sampling
struct EnergyData {
    double zenith;
    double azimuth;
    double meanEnergy;
};

// Global variables to store the muon and energy data
std::vector<MuonData> muon_data;
std::vector<EnergyData> energy_data;

// Function to read the energy data file...i.e. the grid of (zen, az, energy)
void LoadEnergyData(const std::string& filename) {
    energy_data.clear();

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        G4cerr << "Error opening energy data file!" << G4endl;
        return;
    }

    std::string line;
    // Skip the header line
    for (int i = 0; i < 1; ++i) {
        std::getline(infile, line);
    }

    double zenith, azimuth, energy, stddev;
    while (infile >> zenith >> azimuth >> energy) {
        EnergyData data{};
        data.zenith = zenith;
        data.azimuth = azimuth;
        data.meanEnergy = energy;
        energy_data.push_back(data); // adds a new element to the end of energy_data
    }
    infile.close();
}

// Function to find the closest matching energy for given angles
static double FindClosestEnergy(double zenith, double azimuth) {
    double minDistance = 1e10;
    double closestEnergy = 98.0; // Default energy just in case we need it (mean underground energy)

    // Normalize azimuth to [0, 360)...this shouldn't be an issue but it will handle any case where az doesnt belong to [0, 360]
    while (azimuth >= 360.0) azimuth -= 360.0;
    while (azimuth < 0.0) azimuth += 360.0;

    for (const auto& data : energy_data) {
        // Calculate angular distance using both zenith and azimuth
        double dZenith = zenith - data.zenith;
        double dAzimuth = std::min(std::abs(azimuth - data.azimuth),
            360.0 - std::abs(azimuth - data.azimuth));

        // Use Euclidean distance in angle space
        double distance = std::sqrt(dZenith * dZenith + dAzimuth * dAzimuth);

        if (distance < minDistance) {
            minDistance = distance;
            closestEnergy = data.meanEnergy;
        }
    }

    return closestEnergy * GeV; // converting the GeV energies to MeV for geant4......note that the energies in the grid are GeV so we need to do this
}

static void LoadMuonDataAndComputeCDF(const std::string& filename) {

    muon_data.clear(); // cleans up data from previous event to save memory

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        G4cerr << "Error opening data file!" << G4endl;
        return;
    }

    double zenith, azimuth, intensity;
    double total_intensity = 0.0;

    // Read the data from the file
    while (infile >> zenith >> azimuth >> intensity) {
        MuonData data{};
        data.zenith = zenith;
        data.azimuth = azimuth;
        data.intensity = intensity;
        muon_data.push_back(data);
        total_intensity += intensity;  // Accumulate total intensity
    }
    infile.close();

    // Normalize intensities to compute the CDF for sampling
    double cumulative = 0.0;
    for (auto& data : muon_data) {
        data.intensity /= total_intensity;
        cumulative += data.intensity;
        data.cdf = cumulative;
    }
}

// Function to sample from the CDF
static std::pair<double, double> SampleAnglesFromCDF() {
    double random_value = G4UniformRand();

    // Binary search for position in CDF...i.e. is the sampled random number left or right of our mid point
    size_t left = 0;
    size_t right = muon_data.size() - 1;

    // this part checks if left or right and then performs the search based on the result
    while (left < right) {
        size_t mid = (left + right) / 2;
        if (muon_data[mid].cdf < random_value) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }

    // Interpolate between points
    if (left == 0) {
        return std::make_pair(muon_data[0].zenith, muon_data[0].azimuth);
    }

    size_t idx1 = left - 1;
    size_t idx2 = left;

    // Calculate interpolation factor
    double cdf1 = muon_data[idx1].cdf;
    double cdf2 = muon_data[idx2].cdf;
    double alpha = (random_value - cdf1) / (cdf2 - cdf1);

    // Interpolate angles...the initial resolution is very fine, but this is a general approach for more coarse angular grids
    double zenith = muon_data[idx1].zenith + alpha * (muon_data[idx2].zenith - muon_data[idx1].zenith);

    // Special handling for azimuth interpolation to handle 0/360 transition
    double azimuth1 = muon_data[idx1].azimuth;
    double azimuth2 = muon_data[idx2].azimuth;

    // Adjust for circular interpolation if the difference is more than 180 degrees
    if (std::abs(azimuth2 - azimuth1) > 180.0) {
        if (azimuth2 > azimuth1) {
            azimuth1 += 360.0;
        }
        else {
            azimuth2 += 360.0;
        }
    }

    double azimuth = azimuth1 + alpha * (azimuth2 - azimuth1);

    // Normalize azimuth back to [0, 360)
    while (azimuth >= 360.0) azimuth -= 360.0;
    while (azimuth < 0.0) azimuth += 360.0;

    return std::make_pair(zenith, azimuth); // returns the interpolated pair of angles
}

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction() {
    CLHEP::HepRandom::setTheSeed(time(NULL));

    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    // assign muon as primary
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
    fParticleGun->SetParticleDefinition(particle);


    // Load both muon angular distribution and energy data

    // for running from the project root 'muon-induced-secondaries'
    LoadMuonDataAndComputeCDF("./energy_sampling/Mountain_underground_intensities_dd.txt");
    LoadEnergyData("./energy_sampling/Mountain_underground_mean_energies.txt");

    // for running from the release directory with the .exe
    //LoadMuonDataAndComputeCDF("../../energy_sampling/Mountain_underground_intensities_dd.txt");
    //LoadEnergyData("../../energy_sampling/Mountain_underground_mean_energies.txt");
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    if (fTotalAngles.is_open()) {
        fTotalAngles.close();
    }

}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // Sample zenith and azimuthal angles using the CDF
    std::pair<double, double> sampled_angles = SampleAnglesFromCDF();
    double sampled_zenith_deg = sampled_angles.first;    // Keep in degrees for energy lookup
    double sampled_azimuth_deg = sampled_angles.second;  // Keep in degrees for energy lookup

    // Find the corresponding energy for these angles
    double particle_energy = FindClosestEnergy(sampled_zenith_deg, sampled_azimuth_deg); // would comment this out if wanting to set fixed energy
    fParticleGun->SetParticleEnergy(particle_energy); // sets the sampled energy but particle_energy can be replaced with a fixed value if needed...i.e 500 * GeV

    // Convert angles to radians for direction calculation...G4 default units are rad
    double sampled_zenith = sampled_zenith_deg * CLHEP::pi / 180.0;
    double sampled_azimuth = sampled_azimuth_deg * CLHEP::pi / 180.0;

    //defines the radius at which the particles start...just above rock hemisphere
    const G4double sphereRadius = 8.305 * m;

    // samples random theta and phi points on the upper hemisphere
    G4double cosTheta = G4UniformRand();
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    G4double phi = 2.0 * CLHEP::pi * G4UniformRand();

    //computes a random start point on the sphere from the theta and phi values
    G4double x = sphereRadius * sinTheta * std::cos(phi);
    G4double y = sphereRadius * sinTheta * std::sin(phi);
    G4double z = sphereRadius * cosTheta;

    //assigns the position of the particle to particle gun
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

    // assign the sampled angles from CDF
    double cosTheta2 = std::cos(sampled_zenith);
    double sinTheta2 = std::sin(sampled_zenith);
    double cosPhi2 = std::cos(sampled_azimuth);
    double sinPhi2 = std::sin(sampled_azimuth);

    // compute sampled momentum direction from CDF angles...this is where the custom distribution gets used
    G4ThreeVector direction(sinTheta2 * cosPhi2, sinTheta2 * sinPhi2, -cosTheta2);
    fParticleGun->SetParticleMomentumDirection(direction);

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
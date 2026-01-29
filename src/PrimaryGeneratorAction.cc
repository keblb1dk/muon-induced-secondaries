#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
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
#include <algorithm>

//create variables for the angular sampling
struct MuonData {
    double zenith;
    double azimuth;
    double intensity;
    double cdf;
};

//create variables for the mean energy sampling
struct EnergyData {
    double zenith;
    double azimuth;
    double meanEnergy;
};

// Structure for CDF-based energy sampling
struct EnergyCDFData {
    double zenith;
    double azimuth;
    std::vector<double> energies; // GeV
    std::vector<double> cdf;      // [0,1]
};

// Global variables to store the muon and energy data
std::vector<MuonData> muon_data;
std::vector<EnergyData> energy_data;

// Global container for energy CDFs
// Keyed by (zenith, azimuth) in degrees
std::map<std::pair<int, int>, EnergyCDFData> energy_cdf_map;

// Function to load energy CDFs from file
void LoadEnergyCDFs(const std::string& filename) {

    energy_cdf_map.clear();

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        G4cerr << "Error opening energy CDF file: " << filename << G4endl;
        return;
    }

    std::string line;
    int zenith_key = 0;
    int azimuth_key = 0;

    while (std::getline(infile, line)) {

        if (line.empty() || line.find("Energy") != std::string::npos || line[0] == '#') {
            continue;
        }

        if (line.find("Zenith=") != std::string::npos) {

            // Example: Zenith=34.0000, Azimuth=265.0000
            size_t zpos = line.find("Zenith=");
            size_t apos = line.find("Azimuth=");

            double zen = std::stod(line.substr(zpos + 7));

            // Handle azimuth which might be "None" or a number
            std::string az_str = line.substr(apos + 8);
            double az = (az_str.find("None") != std::string::npos) ? -1.0 : std::stod(az_str);

            // Round to nearest integer degree for map key
            zenith_key = static_cast<int>(std::round(zen));
            azimuth_key = (az < 0) ? -1 : static_cast<int>(std::round(az));

            energy_cdf_map[{zenith_key, azimuth_key}] = EnergyCDFData();
            energy_cdf_map[{zenith_key, azimuth_key}].zenith = zen;
            energy_cdf_map[{zenith_key, azimuth_key}].azimuth = az;
            continue;
        }

        // Energy and CDF value
        std::istringstream iss(line);
        double energy_GeV, cdf_val;
        if (!(iss >> energy_GeV >> cdf_val)) {
            continue;
        }

        energy_cdf_map[{zenith_key, azimuth_key}].energies.push_back(energy_GeV);
        energy_cdf_map[{zenith_key, azimuth_key}].cdf.push_back(cdf_val);
    }

    infile.close();

    G4cout << "Loaded energy CDFs for "
        << energy_cdf_map.size()
        << " angular bins." << G4endl;
}

// Function to read the mean energy data file
void LoadEnergyData(const std::string& filename) {
    energy_data.clear();

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        G4cerr << "Error opening energy data file: " << filename << G4endl;
        return;
    }

    std::string line;
    // Skip the header line
    std::getline(infile, line);

    double zenith, azimuth, energy;
    std::string az_str;

    while (infile >> zenith >> az_str >> energy) {
        EnergyData data;
        data.zenith = zenith;
        // Handle "None" azimuth values
        data.azimuth = (az_str == "None") ? -1.0 : std::stod(az_str);
        data.meanEnergy = energy;
        energy_data.push_back(data);
    }
    infile.close();

    G4cout << "Loaded mean energy data for "
        << energy_data.size()
        << " angular bins." << G4endl;
}

// Function to find the closest matching mean energy for given angles
double FindClosestMeanEnergy(double zenith, double azimuth) {
    double minDistance = 1e10;
    double closestEnergy = 288.0; // Default energy

    // Normalize azimuth to [0, 360)
    while (azimuth >= 360.0) azimuth -= 360.0;
    while (azimuth < 0.0) azimuth += 360.0;

    for (const auto& data : energy_data) {
        // Calculate angular distance using both zenith and azimuth
        double dZenith = zenith - data.zenith;

        double dAzimuth;
        if (data.azimuth < 0) {
            // Handle "None" azimuth case (azimuth-independent)
            dAzimuth = 0.0;
        }
        else {
            dAzimuth = std::min(std::abs(azimuth - data.azimuth),
                360.0 - std::abs(azimuth - data.azimuth));
        }

        // Use Euclidean distance in angle space
        double distance = std::sqrt(dZenith * dZenith + dAzimuth * dAzimuth);

        if (distance < minDistance) {
            minDistance = distance;
            closestEnergy = data.meanEnergy;
        }
    }

    return closestEnergy * GeV;
}

// Function to sample energy from CDF for given angles
double SampleEnergyFromCDF(double zenith, double azimuth) {
    // Round angles to nearest integer for map lookup
    int zenith_key = static_cast<int>(std::round(zenith));
    int azimuth_key = (azimuth < 0) ? -1 : static_cast<int>(std::round(azimuth));

    // Normalize azimuth
    while (azimuth_key >= 360) azimuth_key -= 360;
    while (azimuth_key < 0 && azimuth_key != -1) azimuth_key += 360;

    // Try to find exact match
    auto it = energy_cdf_map.find({ zenith_key, azimuth_key });

    // If not found, search for closest bin
    if (it == energy_cdf_map.end()) {
        double minDistance = 1e10;
        std::pair<int, int> closestKey;

        for (const auto& pair : energy_cdf_map) {
            double dZen = zenith - pair.second.zenith;
            double dAz = 0.0;

            if (pair.first.second >= 0 && azimuth_key >= 0) {
                dAz = std::min(std::abs(azimuth_key - pair.first.second),
                    360 - std::abs(azimuth_key - pair.first.second));
            }

            double distance = std::sqrt(dZen * dZen + dAz * dAz);

            if (distance < minDistance) {
                minDistance = distance;
                closestKey = pair.first;
            }
        }

        it = energy_cdf_map.find(closestKey);

        if (it == energy_cdf_map.end()) {
            G4cerr << "Warning: No CDF found for zenith=" << zenith
                << ", azimuth=" << azimuth << ". Using default energy." << G4endl;
            return 288.0 * GeV;
        }
    }

    const EnergyCDFData& cdf_data = it->second;

    // Sample a random number
    double random_value = G4UniformRand();

    // Binary search for the CDF value
    auto lower = std::lower_bound(cdf_data.cdf.begin(), cdf_data.cdf.end(), random_value);

    if (lower == cdf_data.cdf.begin()) {
        return cdf_data.energies[0] * GeV;
    }

    if (lower == cdf_data.cdf.end()) {
        return cdf_data.energies.back() * GeV;
    }

    // Linear interpolation between CDF points
    size_t idx = std::distance(cdf_data.cdf.begin(), lower);
    size_t idx1 = idx - 1;
    size_t idx2 = idx;

    double cdf1 = cdf_data.cdf[idx1];
    double cdf2 = cdf_data.cdf[idx2];
    double E1 = cdf_data.energies[idx1];
    double E2 = cdf_data.energies[idx2];

    // Interpolation factor
    double alpha = (random_value - cdf1) / (cdf2 - cdf1);

    // Interpolate in log space for better accuracy across orders of magnitude
    double logE = std::log(E1) + alpha * (std::log(E2) - std::log(E1));

    return std::exp(logE) * GeV;
}

void LoadMuonDataAndComputeCDF(const std::string& filename) {

    muon_data.clear();

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        G4cerr << "Error opening muon angular data file: " << filename << G4endl;
        return;
    }

    double zenith, azimuth, intensity;
    double total_intensity = 0.0;

    // Read the data from the file
    while (infile >> zenith >> azimuth >> intensity) {
        MuonData data;
        data.zenith = zenith;
        data.azimuth = azimuth;
        data.intensity = intensity;
        muon_data.push_back(data);
        total_intensity += intensity;
    }
    infile.close();

    // Normalize intensities to compute the CDF for sampling
    double cumulative = 0.0;
    for (auto& data : muon_data) {
        data.intensity /= total_intensity;
        cumulative += data.intensity;
        data.cdf = cumulative;
    }

    G4cout << "Loaded angular distribution with "
        << muon_data.size()
        << " angular bins." << G4endl;
}

// Function to sample from the angular CDF
std::pair<double, double> SampleAnglesFromCDF() {
    double random_value = G4UniformRand();

    // Binary search for position in CDF
    size_t left = 0;
    size_t right = muon_data.size() - 1;

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

    // Interpolate angles
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

    return std::make_pair(zenith, azimuth);
}


//G4String samplingDir = "./energy sampling/";
G4String samplingDir = "../../energy sampling/";
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
    fParticleGun(nullptr),
    fMessenger(nullptr),
    fUseMeanSampling(true),  // Default to mean energy sampling
    fEnergyThreshold(0.0),
    fAngularDataFile(""),
    fMeanEnergyFile(""),
    fCDFEnergyFile("")
{
    G4cout << "PrimaryGeneratorAction: Energy threshold initialized to "
        << fEnergyThreshold / GeV << " GeV" << G4endl;

    CLHEP::HepRandom::setTheSeed(time(NULL));

    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
    fParticleGun->SetParticleDefinition(particle);

    // Create messenger for macro control
    fMessenger = new PrimaryGeneratorMessenger(this);

    // Set default file paths (update these to your actual paths)
    //default files to load on launch
    fAngularDataFile = "Test_Mountain_underground_intensities_dd.txt";
    fMeanEnergyFile = "Test_Mountain_underground_mean_energies.txt";
    fCDFEnergyFile = "Test_Mountain_underground_energy_CDFs.txt";

    // Load angular distribution (always needed)
    //LoadMuonDataAndComputeCDF(samplingDir + fAngularDataFile);

    // Load mean energy data by default
    //LoadEnergyData(samplingDir + fMeanEnergyFile);

    // for running from the project root 'muon-induced-secondaries'
    LoadMuonDataAndComputeCDF("../../energy_sampling/Test_Mountain_underground_intensities_dd.txt");
    LoadEnergyData("../../energy_sampling/Test_Mountain_underground_mean_energies.txt");

    // for running from the release directory with the .exe
    //LoadMuonDataAndComputeCDF("./energy_sampling/Test_Mountain_underground_intensities_dd.txt");
    //LoadEnergyData("./energy_sampling/Test_Mountain_underground_mean_energies.txt");

    G4cout << "==================================================" << G4endl;
    G4cout << "PrimaryGeneratorAction initialized" << G4endl;
    G4cout << "Energy sampling mode: " << (fUseMeanSampling ? "MEAN" : "CDF") << G4endl;
    G4cout << "==================================================" << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fMessenger;
    delete fParticleGun;
    if (fTotalAngles.is_open()) {
        fTotalAngles.close();
    }
}

void PrimaryGeneratorAction::SetUseMeanSampling(G4bool val) {
    fUseMeanSampling = val;

    G4cout << "==================================================" << G4endl;
    G4cout << "Switching to " << (val ? "MEAN" : "CDF") << " energy sampling" << G4endl;

    // Load the appropriate data file if switching modes
    if (val && energy_data.empty()) {
        G4cout << "Loading mean energy data..." << G4endl;
        LoadEnergyData(samplingDir + fMeanEnergyFile);
    }
    else if (!val && energy_cdf_map.empty()) {
        G4cout << "Loading CDF energy data..." << G4endl;
        LoadEnergyCDFs(samplingDir + fCDFEnergyFile);
    }

    G4cout << "==================================================" << G4endl;
}

void PrimaryGeneratorAction::SetAngularDataFile(const G4String& filename) {
    fAngularDataFile = filename;
    LoadMuonDataAndComputeCDF(samplingDir + fAngularDataFile);
}

void PrimaryGeneratorAction::SetMeanEnergyFile(const G4String& filename) {
    fMeanEnergyFile = filename;
    if (fUseMeanSampling) {
        LoadEnergyData(samplingDir + fMeanEnergyFile);
    }
}

void PrimaryGeneratorAction::SetCDFEnergyFile(const G4String& filename) {
    fCDFEnergyFile = filename;
    if (!fUseMeanSampling) {
        LoadEnergyCDFs(samplingDir + fCDFEnergyFile);
    }
}

void PrimaryGeneratorAction::SetMuonEnergyThreshold(G4double eThreshold) {
    fEnergyThreshold = eThreshold;

    G4cout << "==================================================" << G4endl;
    G4cout << "PrimaryGeneratorAction: Energy threshold set to "
        << fEnergyThreshold / GeV << " GeV" << G4endl;

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // Sample zenith and azimuthal angles using the CDF
    std::pair<double, double> sampled_angles = SampleAnglesFromCDF();
    double sampled_zenith_deg = sampled_angles.first;
    double sampled_azimuth_deg = sampled_angles.second;

    // Sample energy based on the selected method
    double particle_energy;

    if (fUseMeanSampling) {
        // Use mean energy sampling (faster)
        particle_energy = FindClosestMeanEnergy(sampled_zenith_deg, sampled_azimuth_deg);
        
    }
    else {
        // Use CDF-based energy sampling (more accurate)
        particle_energy = SampleEnergyFromCDF(sampled_zenith_deg, sampled_azimuth_deg);
    }

    if (particle_energy >= fEnergyThreshold) {
        fParticleGun->SetParticleEnergy(particle_energy);
        //G4cout << "Energy: " << particle_energy / GeV << " GeV" << G4endl;
    }
    // Convert angles to radians for direction calculation
    double sampled_zenith = sampled_zenith_deg * CLHEP::pi / 180.0;
    double sampled_azimuth = sampled_azimuth_deg * CLHEP::pi / 180.0;

    // Define the radius at which the particles start
    const G4double sphereRadius = 8.305 * m;

    // Sample random theta and phi points on the upper hemisphere
    G4double cosTheta = G4UniformRand();
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    G4double phi = 2.0 * CLHEP::pi * G4UniformRand();

    // Compute a random start point on the sphere from the theta and phi values
    G4double x = sphereRadius * sinTheta * std::cos(phi);
    G4double y = sphereRadius * sinTheta * std::sin(phi);
    G4double z = sphereRadius * cosTheta;

    // Assign the position of the particle to particle gun
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

    // Assign the sampled angles from CDF
    double cosTheta2 = std::cos(sampled_zenith);
    double sinTheta2 = std::sin(sampled_zenith);
    double cosPhi2 = std::cos(sampled_azimuth);
    double sinPhi2 = std::sin(sampled_azimuth);

    // Compute sampled momentum direction from CDF angles
    G4ThreeVector direction(sinTheta2 * cosPhi2, sinTheta2 * sinPhi2, -cosTheta2);
    fParticleGun->SetParticleMomentumDirection(direction);

    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Write the angular and energy sampling results to a file (if enabled)
    if (fTotalAngles.is_open()) {
        fTotalAngles << sampled_zenith * 180 / CLHEP::pi << " "
            << sampled_azimuth * 180 / CLHEP::pi << " "
            << particle_energy / GeV << std::endl;
    }
}

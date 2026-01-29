# Muon-Induced Secondary Particle Simulation

[![DOI](https://doi.org/10.5281/zenodo.18408949.svg)](https://doi.org/10.5281/zenodo.18408949)


This project uses [Geant4](https://geant4.web.cern.ch/) to simulate muon-induced secondary particles in underground laboratories. It combines realistic muon flux distributions from the [MUTE package](https://github.com/wjwoodley/mute) with Monte Carlo particle tracking to study secondary particle production in underground detector environments.

## Project Overview

The simulation models cosmic muons penetrating an underground laboratory situated within a rock hemisphere. The lab geometry consists of a 3.05×5.2×3.05 meter cavity. Muons are sampled from realistic angular and energy distributions, and secondary particles created through muon interactions in the surrounding rock are tracked as they enter the laboratory space.

### Key Features

- Realistic muon flux sampling using MUTE-generated distributions
- Underground laboratory geometry in rock hemisphere  
- Secondary particle tracking and analysis
- Angle and energy-dependent muon sampling
- Detailed particle interaction data output

## Project Structure

```
muon-induced-secondaries/
├
├── src/                            # Source code files
├── include/                        # Header files  
├── energy_sampling/                # Input data files
│   ├── Mountain_underground_mean_energies.txt
│   ├── Mountain_underground_intensities_dd.txt  
│   ├── Mountain_underground_flux_grid.txt
|   ├── example_mountain_profile.txt
│   └── ....
|   
├── macros/                         # Geant4 macro files
├── secondaries/                    # Output directory
│   └── secondary_data/
│       └── escaped_secondaries_[LabName].csv
├── Rock Composition Data/                         # contains chemical rock composition txt files for different labs
├── main.cc                         # Main program entry point
├── CMakeLists.txt                  # Build configuration
├── calculate_e_sampling.py         # Python preprocessing script
├── mute_update.txt                 # function for MUTE to generate angle-dependent energy grid
└── README.md                       # This file
```

## Prerequisites

- **Geant4** (version 10.7 or higher recommended)
- **CMake** (version 3.5 or higher)
- **C++ compiler** with C++11 support
- **Python 3** with NumPy and SciPy
- **MUTE package** for muon flux calculations

## Installation and Setup

### 1. Install MUTE Package

**Recommended: Install Modified MUTE Fork**

Install the modified MUTE package directly from the underground-energy-sampling branch:
```bash
pip install git+https://github.com/keblb1dk/mute.git@underground-energy-sampling
```

This fork includes the necessary modifications for flux grid output and angle-dependent energy sampling.
For this option, if you already have MUTE installed, please first do
```bash
pip uninstall mute
```

**2. Alternative: Install Standard MUTE and Modify Manually**

If you prefer to use the standard MUTE package and apply modifications yourself:

1. Install the standard MUTE package following the instructions at: https://github.com/wjwoodley/mute

2. Modify the MUTE source code to enable flux grid output:
   - Navigate to your MUTE installation directory (typically located at):
```
     Users/[username]/AppData/Local/Programs/Python/Python312/Lib/site-packages/mute/underground.py
```
   - Open `underground.py` and find the `calc_u_e_spect()` function
   - Copy the code from `mute_update.txt` in this project and completely replace the calc_u_e_spect() function and add the two new functions, calc_u_e_cdf_sampling() and calc_u_e_sampling, directly below it.

Note: The manual modification approach requires updating the code each time MUTE is reinstalled or updated.

### 3. Build the Geant4 Simulation

Navigate to the project root directory and build:

```bash
cd muon-induced-secondaries
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

This creates the `MuSecondaries` executable in the build directory.

## Usage

### Step 1: Generate Input Data (Optional)

The project comes with pre-generated example data from a mountain profile. To generate new data for a different location or configuration:
```bash
cd muon-induced-secondaries
python calculate_e_sampling.py
```
The script runs in three phases:
Phase 1: Generate MUTE Data

The script will use MUTE to calculate underground muon flux distributions and create two key files:

- [LabName]_underground_flux_grid.txt - Full energy-angle flux distributions
- [LabName]_underground_intensities_dd.txt - Angular intensity weights

These files are saved to your MUTE data directory (typically):

```Users/[username]/AppData/Local/Programs/Python/Python312/Lib/site-packages/mute/data/underground/```

Important: The script will pause here and wait for you to press Enter. During this pause, you must:

- Navigate to the MUTE data directory above
- Copy the two generated files
- Paste them into your muon-induced-secondaries/energy_sampling/ directory
- Return to the terminal and press Enter to continue

Phase 2: Generate Mean Energies

After copying the files, the script will:

Process the flux grid to calculate angle-dependent mean energies or angle-dependent CDFs for energy sampling

Create ```[LabName]_underground_mean_energies.txt``` or ```[LabName]_underground_energy_CDFs.txt```in the energy_sampling/ directory

Phase 3: Validation (Optional)

The script will pause again, asking if you want to compare global mean energies calculated using two different integration methods. This validates the accuracy of the angle-dependent approach against the direct integrated spectrum method.

Note: The lab name is set by ```mtc.set_lab("Mountain")``` in the script. Change this to match your desired configuration, and all generated files will use your specified lab name.

### Step 2: Run the Simulation

From the project root directory:

```bash
# Interactive mode (with visualization)
./build/Release/MuSecondaries

# Batch mode with macro file
./build/Release/MuSecondaries macros/your_macro.mac
```

**Important**: Always run from the `muon-induced-secondaries` root directory, as the simulation uses relative paths to locate input files.

### Step 3: Analyze Output

Results are saved to `secondaries/secondary_data/escaped_secondaries.csv` and include:

- Track and event IDs
- Particle creation and exit energies
- Position coordinates (creation and exit points)
- Momentum directions
- Rock depth traversed
- Primary muon properties
- Creation process information

## Customization

### Laboratory Dimensions

The laboratory cavity and rock hemisphere dimensions are defined in two files and can be updated manually:

**DetectorConstruction.cc**:
```cpp
G4double labX = 3.05 / 2. * m;  // Half-width in X
G4double labY = 5.2 / 2. * m;  // Half-width in Y  
G4double labZ = 3.05 / 2. * m;  // Half-height in Z
```

**SteppingAction.cc**:
```cpp
G4double labX = 3.05 / 2.0;     // Half-width in X (meters)
G4double labY = 5.2 / 2.0;     // Half-width in Y (meters)
G4double labZ = 3.05;           // Height in Z (meters)
```

### Rock and World Geometry

Hemisphere dimensions in **DetectorConstruction.cc**:
```cpp
G4double outerRadius = 8.4 * m;      // World volume radius
G4double outerRadiusRock = 8.3 * m;  // Rock hemisphere radius
```

### Input Data Files

Replace files in `energy_sampling/` with your own data:
- `Mountain_underground_mean_energies.txt` - Angle-dependent mean energies
- `Mountain_underground_energy_CDFs.txt` - Angle-dependent energy CDFs for more accurate sampling
- `Mountain_underground_intensities_dd.txt` - Angular intensity distribution
- `Mountain_underground_flux_grid.txt` - Full energy-angle flux grid
- `example_mountain_profile.txt` - mountain profile grid required by MUTE: (Theta, Azimuth, Slant Depth)

When using custom data files, ensure you update the `PrimaryGeneratorAction.cc` file with the correct intensity and energy file names. Please look at the included example files to ensure your formatting is correct.
More detailed information on how set files and other observables can be found in [ReadMe](https://github.com/keblb1dk/muon-induced-secondaries/blob/main/src/ReadMe.md).

## Troubleshooting

### Common Issues

**Build errors**: Ensure Geant4 environment is properly sourced and CMake can find Geant4 installation

**File not found errors**: Verify you're running from the `muon-induced-secondaries` root directory

**Empty output**: Check that input files exist in `energy_sampling/` directory

**Python import errors**: Ensure MUTE package is properly installed and the source modification was applied

### Support

For MUTE-related issues, consult the official documentation: https://github.com/wjwoodley/mute

For Geant4 support, refer to the Geant4 User's Guide: https://geant4-userdoc.web.cern.ch/

## Citation
The current version can be cited with the [Zenodo DOI](https://doi.org/10.5281/zenodo.18408949).
Please also cite the APS PRD [article](https://doi.org/10.1103/7n8l-pr2g) for this work.

The citations for MUTE and Geant4 can be found at [MUTE](https://github.com/wjwoodley/mute) and [Geant4](https://geant4.web.cern.ch/).

## Authors
[Dakota Keblbeck](dakota_keblbeck@mines.edu)

## License

This project builds upon the MUTE package and Geant4 toolkit. Please respect the licenses of these dependencies when using or redistributing this code.

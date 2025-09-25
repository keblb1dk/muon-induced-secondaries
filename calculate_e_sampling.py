import mute.constants as mtc
import mute.underground as mtu
import numpy as np




#Example of how to use MUTE to calculate the energy spectrum and mean energy, and the underground intensities
#Please see  https://github.com/wjwoodley/mute/tree/main  for in-depth examples

mtc.set_verbose(2)
mtc.set_lab("Mountain")
lab_name = mtc.get_lab()
mtc.set_overburden("mountain")
mtc.set_medium("rock")
mtc.set_reference_density(2.65)
mtc.set_n_muon(1000000)

# replace with custom file as needed
mtc.load_mountain("./energy_sampling/example_mountain_profile.txt")


#The energy spectrum, flux grid, and intensities will usually get saved to  Users\myUserName\AppData\Local\Programs\Python\Python312\Lib\site-packages\mute\data\underground
# these files need to be copied into the directory muon-induced-secondaries/energy_sampling/

mountain_u_e_spect = mtu.calc_u_e_spect(output=True, model = "daemonflux")  # calculates energy spectrum and saves the flux grid
mtu.calc_u_intensities(method = "dd", output = True, model = "daemonflux")  # calculates underground intensities

mountain_u_energies = mtu.calc_u_mean_e(model = "daemonflux") # calculate mean energy
mean_energy = mountain_u_energies[0]

print("The next section takes the flux grid and computes angle-dependent mean energies")
input("Press Enter to continue...")


def compute_mean_energies(input_file, output_file):
    results = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    current_zenith = None
    current_azimuth = None
    energies = []
    intensities = []
    
    for line in lines:
        line = line.strip()
        
        if not line or line.startswith('Energy'):
            continue
            
        if line.startswith('Zenith='):
            # Process previous section
            if energies and np.sum(intensities) > 0:
                energies_arr = np.array(energies)
                intensities_arr = np.array(intensities)
                norm_intensities = intensities_arr / np.trapz(intensities_arr, energies_arr)
                mean_energy_mev = np.trapz(energies_arr * norm_intensities, energies_arr)
                mean_energy_gev = mean_energy_mev / 1000.0
                results.append((current_zenith, current_azimuth, mean_energy_gev))
            
            # Parse new angles
            parts = line.split(', ')
            current_zenith = float(parts[0].split('=')[1])
            current_azimuth = float(parts[1].split('=')[1])
            energies = []
            intensities = []
            
        else:
            # Parse data
            try:
                energy, intensity = line.split()
                energies.append(float(energy))
                intensities.append(float(intensity))
            except ValueError:
                continue
    
    # Process final section
    if energies and np.sum(intensities) > 0:
        energies_arr = np.array(energies)
        intensities_arr = np.array(intensities)
        norm_intensities = intensities_arr / np.trapz(intensities_arr, energies_arr)
        mean_energy_mev = np.trapz(energies_arr * norm_intensities, energies_arr)
        mean_energy_gev = mean_energy_mev / 1000.0
        results.append((current_zenith, current_azimuth, mean_energy_gev))
    
    # Write output
    with open(output_file, 'w') as f:
        f.write("Zenith(degrees)\tAzimuth(degrees)\tMeanEnergy(GeV)\n")
        for zenith, azimuth, mean_energy in results:
            f.write(f"{zenith:.4f}\t{azimuth:.4f}\t{mean_energy:.4f}\n")
    
    print(f"Processed {len(results)} angle pairs")


input_file = f"./energy_sampling/{lab_name}_underground_flux_grid.txt"
output_file = f"./energy_sampling/{lab_name}_underground_mean_energies.txt"
compute_mean_energies(input_file, output_file)


print("The next section compares the mean energy from the integrated spectrum to that calculated using the mean energy grid")
input("Press Enter to continue...")


import numpy as np
from scipy import integrate as scii

def compute_global_mean_intensity_weighted(mean_energy_file, flux_grid_file):
    """
    Compute global mean energy weighted by the total intensity at each angle
    Uses direct integration without interpolation
    """
    # Read the mean energy data
    mean_data = np.loadtxt(mean_energy_file, skiprows=1)
    zenith_deg = mean_data[:, 0]
    azimuth_deg = mean_data[:, 1] 
    mean_energies = mean_data[:, 2]  # Already in GeV
    
    # Read flux grid to get intensities using integration over energy
    angle_intensities = {}  # Dictionary to store (zenith, azimuth): total_intensity
    current_zenith = None
    current_azimuth = None
    current_energies = []
    current_intensities = []
    
    with open(flux_grid_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('Energy'):
                continue
                
            if line.startswith('Zenith='):
                # Save previous total using integration if exists
                if current_zenith is not None and current_energies and current_intensities:
                    total_flux = scii.simpson(y = current_intensities,x =  current_energies)
                    angle_intensities[(current_zenith, current_azimuth)] = total_flux
                
                # Parse new angles
                parts = line.split(', ')
                current_zenith = float(parts[0].split('=')[1])
                current_azimuth = float(parts[1].split('=')[1])
                current_energies = []
                current_intensities = []
                
            else:
                # Collect energy and intensity data for integration
                try:
                    energy, intensity = line.split()
                    current_energies.append(float(energy))
                    current_intensities.append(float(intensity))
                except ValueError:
                    continue
        
        
        if current_zenith is not None and current_energies and current_intensities:
            total_flux = scii.simpson(y = current_intensities,x =  current_energies)
            angle_intensities[(current_zenith, current_azimuth)] = total_flux
    
    # Match intensities with mean energies (only use points that exist in both)
    matched_zen = []
    matched_az = []
    matched_energies = []
    matched_intensities = []
    
    for i in range(len(zenith_deg)):
        key = (zenith_deg[i], azimuth_deg[i])
        if key in angle_intensities and angle_intensities[key] > 0:
            matched_zen.append(zenith_deg[i])
            matched_az.append(azimuth_deg[i])
            matched_energies.append(mean_energies[i])
            matched_intensities.append(angle_intensities[key])
    
    # Convert to arrays
    matched_zen = np.array(matched_zen)
    matched_az = np.array(matched_az)
    matched_energies = np.array(matched_energies)
    matched_intensities = np.array(matched_intensities)
    
    print(f"Data points matched: {len(matched_energies)}")
    print(f"Zenith range: {np.min(matched_zen):.1f}° to {np.max(matched_zen):.1f}°")
    print(f"Energy range: {np.min(matched_energies):.2f} to {np.max(matched_energies):.2f} GeV")
    
    # Group data by zenith angle for direct integration (no interpolation)
    zenith_groups = {}
    for i in range(len(matched_zen)):
        zen = matched_zen[i]
        if zen not in zenith_groups:
            zenith_groups[zen] = {'az': [], 'energy': [], 'intensity': []}
        zenith_groups[zen]['az'].append(matched_az[i])
        zenith_groups[zen]['energy'].append(matched_energies[i])
        zenith_groups[zen]['intensity'].append(matched_intensities[i])
    
    # Integrate over azimuth for each zenith
    zenith_values = []
    azimuth_integrated_weighted = []
    azimuth_integrated_intensity = []
    
    for zen in sorted(zenith_groups.keys()):
        data = zenith_groups[zen]
        az_sorted_idx = np.argsort(data['az'])
        
        # Sort all arrays by azimuth
        weighted_energy_az = np.array(data['energy'])[az_sorted_idx] * np.array(data['intensity'])[az_sorted_idx]
        intensity_az = np.array(data['intensity'])[az_sorted_idx]
        az_sorted = np.array(data['az'])[az_sorted_idx]
        
        # Integrate over azimuth 
        weighted_integral = scii.simpson(y = weighted_energy_az, x=np.radians(az_sorted))
        intensity_integral = scii.simpson(y = intensity_az, x=np.radians(az_sorted))
        
        zenith_values.append(zen)
        azimuth_integrated_weighted.append(weighted_integral)
        azimuth_integrated_intensity.append(intensity_integral)
    
    # Convert to arrays for zenith integration
    zenith_values = np.array(zenith_values)
    azimuth_integrated_weighted = np.array(azimuth_integrated_weighted)
    azimuth_integrated_intensity = np.array(azimuth_integrated_intensity)
    
    # Integrate over zenith with cos(zenith) weighting 
    # cos(zenith) is decreasing, so take absolute value 
    total_weighted = abs(scii.simpson(y = azimuth_integrated_weighted, x=np.cos(np.radians(zenith_values))))
    total_intensity = abs(scii.simpson(y = azimuth_integrated_intensity, x=np.cos(np.radians(zenith_values))))
    
    global_mean_gev = total_weighted / total_intensity
    
    return global_mean_gev * 1000  # Convert to MeV for comparison


if __name__ == "__main__":
    mean_energy_file = f"./energy_sampling/{lab_name}_underground_mean_energies.txt"
    flux_grid_file = f"./energy_sampling/{lab_name}_underground_flux_grid.txt"
    
    print(f"Processing {mean_energy_file}")
    
    try:
        intensity_weighted_mean = compute_global_mean_intensity_weighted(mean_energy_file, flux_grid_file)
        print(f"Intensity-weighted integration: {intensity_weighted_mean:.2f} MeV")
        print(f"Expected from direct calc_u_mean_e: {mean_energy:.2f} MeV")
        
        difference_percent = abs(mean_energy - intensity_weighted_mean) / mean_energy * 100
        print(f"Difference: {difference_percent:.4f}%")
        
    except Exception as e:
        print(f"Intensity-weighted integration failed: {e}")
        import traceback
        traceback.print_exc()



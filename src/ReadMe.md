# User-defined options in command line

There are three main files which can be used by the user in the command line or through a macro file to set observables.

## SteppingAction
The SteppingAction.cc file accepts [LabName] from DetectorMessenger to update the output csv file. This can be done through the command

```/detector/setLabName YourLabName```

This command **must** come before the ```/run/initialize``` command, and can only be set once.

## DetectorConstruction
The DetectorConstruction allows the user to set the density of the rock and to set the rock chemical composition

- set density (default units are g/cm3)
  ```/detector/setDensity value```
- set composition using files in Rock Composition Data
  ```/detector/setRockComp labname_rock.txt```

## PrimaryGeneratorAction
PrimaryGeneratorAction allows you to set an energy threshold, change from mean to CDF sampling, and set the sampling files without rebuilding the executable.

- set energy threshold (default value and units are 0 and GeV)
  ```/primary/setEnergyThreshold 0.0```

- change from mean energy sampling (default) to CDF sampling
  ```/primary/useMeanSampling false```

  When you change between mean and CDF sampling you **must** update the sampling files accordingly.

- set the angular sampling file
  ```/detector/setAngularFile angular_file.txt```

- set the mean energy sampling file
  ```/primary/setMeanEnergyFile mean_energy_file.txt```

- set the mean energy sampling file
  ```/primary/setCDFEnergyFile mean_energy_file.txt```

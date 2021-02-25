# Geant4 Zero Knowledge Bremstrahhlung Production Simulation

Name: Readme.md  
Date: 25 August 2018  
Auth: J. Vavrek and B.S. Henderson  
Mail: jvavrek@mit.edu  

## Usage

See below for build instructions

ZKBrem (vis option) (macro name) (seed/label) (LYSO/LaBr3) (plate/target) (beam E) (beam dist) (particle) (false dets) (filters) (beamoffset)

- ZKBrem with no options provided will print these usage instructions
- First option (vis option) is mandatory to pass the help message, once vis
  option is set defaults proceed as follows, must supply all preceding arguments
- Defaults: ZKBrem (vis option) ZKBREM_TOPDIR/runtime/ZK.mac 0 0 0 0 0 0 false true 0
- vis option: visOff (no visualization), visOn (default Geant4 visualization),
  visQt (Qt visualization)
- macro name: Full path to Geant4 macro, will not be used in Qt visualization
- Seed/label: Integer used to label output file and offset the random seed
- LYSO/LaBr3: In-beam scintillator/shielding configuration, see documentation
  noted at instances of "lysoconf" in src/geometryConstruction.cc
- plate/target: Foil/target configurations, see documentation noted at instances
  of "plateconf" in src/geometryConstruction.cc, uses bitwise test for the
  Sept. 2017 configurations
- beam E: Beam energy parameter, 0 = Uniform energy between 1.9 and 2.3 MeV
  (range can be changed in src/PGA.cc), 1 = Sampling of U-238/Al-27 resonances
  2 = Sept. 2017 endpoint energy fixed (2.521 MeV)
- beam dist: Beam illumination mode, 0 = collimated bremsstrahlung beam,
  1 = illumination of HPGe detectors originating from the foils (in their Sept.
  2017 configuration), 2 = electron beam on radiator for full brem generation
- particle: Beam particle type, 0 = gamma, 1 = geantino, 2 = electron
- false dets: boolean for inclusion of false air detector planes immediately
  after the radiator and before the foil used for recording the bremsstrahlung
  production, otherwise produces masssive output files
- filters: boolean for inclusion of lead filters between HPGe detectors and the
  foil(s), useful for simulations of intrinsic detector efficiencies
- beamoffset: Parameter passed to the PGA for utility, currently unused

Example configurations:

ZKBrem visQt runtime/ZK.mac 1 4 768 1 1 0 0 1
- "fullsim" from the NRF papers, in which the photons origin from within the
  foils, and detectors are fully simulated

ZKBrem visQt runtime/ZK.mac 1 4 768 1 0 0 0 1
- "brute force" from the NRF papers in which the collimated beam is simualted,
  add appropriate target code to the 768 plate/target argument to add a
  transmission target to the standard foils

ZKBrem visQt runtime/ZK.mac 1 4 4 1 1 0 0 0
- intrinsic efficiency simulation for the HPGes to get right illumination with
  filters and physical foils removed

ZKBrem visQt runtime/ZK.mac 1 4 768 1 1 1 0 0
- geometric efficiency simulation for the HPGes to get right illumination with
  geantinos

ZKBrem visQt runtime/ZK.mac 1 4 768 2 2 2 1 1
- generation of bremsstrahlung from the electron beam (fixed endpoint set in
  source code)



## Overview 

ZKBrem is designed to provide the ability to rapidly generate
bremstrahhlung emission data (energy spectra, angular distributions,
etc) for electron beams into various target materials and store them
in ROOT files for later use as particle sources distributions for
further simulation. The simulation is derived from ZKTemplate, and,
therefore, incorporates several enhanced capabilities including
optional MPI parallelization and ROOT input/output.


## Directory Structure

ZKBrem has the following directory structure:

- **ZKBrem.cc:** standard G4 top-level file containing mandatory C++ main() function.

- **makefile:** GNU makefile for building ZKTemplate

- **Readme.md:** this file

- **include/:** standard class header files (*.hh)

- **src/:** : standard class source files (*.cc)

- **lib/:** contains the ROOT library and dictionary files
 
- **runtime/:** contains G4GPS, visualization, and other macros


## Installation

The user's should first configure their environment and compile ZKBrem
to ensure everythign is working correctly. To do so:

 1. Ensure your system contains an up-to-date, complete installation
 of Geant4 (http://geant4.cern.ch/) with the extended libraries. It is
 preferable (although not optional) to build the Qt libraries and
 include them when building the Geant4 libraries. Ensure that your
 environment and G4 working directory paths are correctly configured.

 2. Ensure your system contains an up-to-date, complete installation
 of the ROOT data analysis toolkit (http://root.cern.ch/drupal/) and
 that your environment has been correctly configured.

 3. If you desire to build simulations with parallel processing,
 ensure your system contains an up-to-date, complete installation of
 MPI ([Open MPI](http://www.open-mpi.org/) or
 [MPICH2](http://www.mcs.anl.gov/research/projects/mpich2staging/goodell/index.php)
 has been used successfully with my custom MPI interface for Geant4.)
 To install MPI easily on MacOS:
 ```bashrc
 $ brew install openmpi
 ```

 4. Unpack the simulation tarball into your geant4 working directory

 5. Add the following lines to your .bashrc file (typically found at
   /home/$USERNAME/.bashrc):
   ```bash
   source /path/to/ZKBrem.setup.sh >& /dev/null
   export ZKBREM_MPIHOME=/absolute/path/to/MPI/installation
   ```

## Building the code

Buildilng the ZKBrem simulation is presently handled by a GNU makefile
located in the top-level directory. To build the sequential binary
called *ZKBrem*:
```bash
$ . ./setup.sh
$ make -jX
```
where 'X' is number of cores to use in the compilation process. 

To build the parallel binary called *ZKBrem_MPI*:
```bash
./parallel.sh
```
This will run a bashrc script that controls a number of internal
variable and environment settings and then calls the GNU makefile with
the correct options. Note that in parallel, the binary name is
"ZKBrem_MPI" Note: this requires that openmpi be installed.

Typing 'make clean' will clean up all the transient build files while
'make libclean' will clean up the ROOT library and dictionary build
files contained in the lib/ directory.


## Executing the code

To run the sequential version from any directory, simply type:
```bash
$ ZKBrem <arg1> <arg2>
```
where arg1 and arg2 are optional cmd line arguments. arg1 can take the
following values: "visOff" (no visualization), "visOn" visualization
with standard OpenGL, "visQt" visualization with enhanced Qt
interface.  arg2 specifies the name of a macro file to use in parallel
processing.

To run the parallel version from any directory, standard MPI methods
must be used. Note that for parallel runs, the "visOff" cmd must be
passed as the first cmd line argument. For example, to run ZKBrem on
12 processors, one would type:
```bash
$ mpirun -np 12 ZKBrem_MPI visOff runtime/ZK.PrototypeRun.mac
```
where "runtime/ZK.PrototypeRun.mac" is an optional argument to specify
a macro file containing a series of G4 commands to run. (If the user
does not specify a parallel macro, the simulation will automatically
run commands contained in the macro "runtime/ZK.mpi.mac".)

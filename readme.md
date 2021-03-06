# fABMACS
A fast, hard-wired version of [GROMACS 5.0.5] that computes free energy estimates via adaptive biasing potentials. The accompanying manuscript can be found at http://scitation.aip.org/content/aip/journal/jcp/145/15/10.1063/1.4964776. Usage tutorials are currently being written to guide new fABMACS users  and can be found at our [wiki](https://github.com/BradleyDickson/fABMACS/wiki).

This code is modified for the following functions:

- **NEW:** Compute escape times using [Hyperdynamics](http://public.lanl.gov/afv/afv-66-hyperPRL.pdf)

- Compute free energy via mABP in 2 RMSD collective variables

- Compute Phi-Psi free energy for ALANINE DIPEPTIDE via mABP and WTmetaD

- Integrative patching script can customize the code to your application

**This is not standard gromacs, don't use it for equilibration tasks!**

All the things required for building the published ligand simulation inputs are included in the [RUNdirs] subdirectory. Alanine dipeptide inputs are included there as well, just check the sub directory names and read the README files in those directories.

The [RUNdirs] directory contains parameter files for everything that appeared in publication. The bias parameter inputs are largely transferable but may need to be tweaked for some situations.

We figure you are in a UNIX environment.

Lastly, distance units are in nanometers, as per GROMACS.

**Rectangular simulation cells are now supported.**
-the "width" variable in params.in now **requires** three (3) "widths" 
-the input now must be "widthx widthy widthz" as reflected in the [RUNdirs] inputs

# To Build:
1. Go to your fABMACS directory. (you've already downloaded, unpacked, etc...)

2. Prepare your "list" file. This file holds the atom number of each atom used in the CV. This release uses RMSD for CV, so no atom is expected to appear twice. The "list" file used for simulations that were published is found in the fABMACS dir and is named "benzoisoxazoloazepinelists". To re-run published simulations, copy this file to the name "list". 

3. Patch the code by running PATCHscript.sh. For Alanine dipeptide ```./PATCHscript.sh alanine``` For custom simulations: Run ./PATCHscript.sh BMAX NPARTS NCV1 NCV2 CVMAX CVREST

 - BMAX is the number of bins that will span the CV space
 - NPARTS is the total number of atoms in the collective variables
 - NCV1 is the number of atoms in CV1 (the first CV)
 - NCV2 is the number of atoms in the second CV
 - CVMAX is the maximum possible value for any CV
 - CVREST is the starting position of a harmonic restraint acting on CV1 and CV2
 - **Example use:** Published ligand simulations used ```./PATCHscript.sh 480 8 4 4 6 5.5``` to run with spherical restraint
 - **Example use:** Published ligand simulations used ```./PATCHscript.sh 480 8 4 4 6 3 cylinder``` to run with cylindrical restraint. CVREST=3 defines the length of the cylinder.
 - **Alanine with Hyperdynamics:** Patch as ```./PATCHscript.sh alanine HYPER``` and see [RUNdirs] for proper params.in files.
 - **Alanine with Overfill limits:** Patch as ```./PATCHscript.sh alanine OVERFILL``` and see [RUNdirs] for proper params.in files.
<a name="patchhyper"></a>
 - **Ligand escape with Hyperdynamics:** Published ligand hyperdynamics simulations used ```./PATCHscript.sh 480 8 4 4 6 3 cylinder HYPER``` where bias, cylinder, and topology inputs are within [RUNdirs]. See [below](#hyperdetail) for further details.

4. Configure GROMACS build as usual. Make your Build and Bin directories. Then run the cmake command in your build dir with the options you need to use to compile standard GROMACS 5.0.5. We used the options ```-DGMX_BUILD_OWN_FFTW=ON  -DGMX_SIMD=AVX2_256 -DGMX_OPENMP=OFF -DGMX_MPI=ON``` Our cmake looked like this:

 - ```cmake PATH-TO-SOURCE -DGMX_BUILD_OWN_FFTW=ON  -DGMX_SIMD=AVX2_256 -DGMX_OPENMP=OFF -DGMX_MPI=ON -DCMAKE_INSTALL_PREFIX=PATH-TO-BIN```

5. Run make from the fABMACS build directory

6. Run make install from the fABMACS build directory


# To run simulations of alanine dipeptide

0. You must have used ```./PATCHscript.sh alanine``` to build fABMACS

1. Go to fABMACS/RUNdirs/ALANINE and build a new tpr file:
 - PATH-TO-grompp_mpi -f md.mdp -c isob.gro -t isob.cpt -o ala.tpr

2. Edit the params.in file to your liking, and make sure you run using the executable that was built using fABMACS. The alanine simulations require the file called "reffreeE" which holds the reference free energy for computing convergence profiles. (see outputs below) "refreeE" is in the RUNdirs/ALANINE directory.

 - Our alanine dipeptide simulations used 8 core and ran as ```mpirun ./mdrun -deffnm``` where we used a symbolic link to define mdrun.

3. You can adjust the parameters and see how things change. 

# Alanine dipeptide Outputs

1. The simulations will write a file named "freeE" that contains the current free energy estimate. The Phi-Psi angles are given in the first two columns, the free energy estimate is given in the third column.

2. Simulations also write a file named "fort.88" The first column is timestep, second and third columns are collective variables (angles), the fourth column is the convergence metric (equation 28 [HERE])

# Custom simulation or re-run our ligand simulations for *free energy*
***Things you need, can all be found in RUNdirs/RErun directory***

- Reference file: Holds position of every atom in the CVs at time t=0. The Reference file for our ligand simulations can be seen in the [RUNdirs] directory named RErun. Your Reference file can be created easily using this bit:
 ```a=`cat PATHto/fABMACS/list`;for w in $a; do grep ' '$w' ' PATHto/EQ.gro |awk '{print $4,$5,$6}';done > Reference```
where "PATHto" is a path to fABMACS where the list file is or a "PATHto" an equilibrated gro file (called EQ.gro here).
- sphpoints file: Holds position of spherical restraint center. The one used in our publication is in the [RUNdirs] directory named RErun. The sphere can be centered anywhere. You need this if you are not using cylindrical restraint. ***Make the radius LARGE if you don't want this restraint to act***
- cylpoints file: Holds two points to define the cylindrical restraint. The one used in our publication is in the [RUNdirs] directory named RErun. We use [VMD] to draw cylinders and select the points. You only need this if you use the cylindrical restraint.
- params.in file: Holds all bias and restraint parameters. This file is described every time the PATCHscript.sh is executed and a general set of parameter values is given. The params.in that were used to run our ligand simulations are in the [RUNdirs]/RErun directory. 

The topolog and equilibrated coordinates (and cpt), and md.mdp are in the RErun directory. Simulation inputs can be built by using: 

***Run from bound state***
 - PATH-TO-grompp_mpi -f md.mdp -c isob.gro -t isob.cpt -o run.tpr

***Run from unbound state (spherical restraint only)***
 - PATH-TO-grompp_mpi -f md.mdp -c short.gro -t short.cpt -o run.tpr

Simulations use RMSD for CVs, so you also need to restrain something in the system so that the reference used to define RMSD is always valid. We do this by adding some restraints via the GROMACS genrestr tool. Be sure to add restraints for you simulations, you can read the topol files for our ligand system to see how we've added these restraints. See the [RUNdirs]/RErun directory, look for back.itp in the topol file.

Running fABMACS simulations is exactly like running standard GROMACS simulations, except that you need the above input files. If you use the cylinder restraint, you need clyploints, otherwise you need sphpoints. Examples of all of these are included, and the params.in file will be suggested and described everytime you patch the code. Logs are also generated when you patch.

Go run simulations! Be sure that you point to the fABMACS executable. Use SPHERE-params.in and CYLINDER-params.in as templates to create your params.in file. If you are re-running our simulations, you can just copy the file that matches the restraint you built in the patching step.

# Simulation Output
1. The simulations will write a file named "freeE" that contains the current free energy estimate. CV1 and CV2 are given in the first two columns, the free energy estimate is given in the third column and the raw sampling histogram is given in the fourth column.

2. Simulations also write a file named "fort.88" The first column is timestep, second and third columns are collective variables 1 and 2, the fourth column is the "hill height"

3. An xyz file of the atoms in the CVs is output, named "fort.81." This file can be used to check that the periodic boundaries are treated correctly.


# <a name="hyperdetail"></a> Custom simulation or re-run our ligand simulations for *Hyperdynamics*
***Things you need, can all be found in RUNdirs/RErun directory***

*The patching command for hyperdynamics simulations was given [above](#patchhyper).*

When hyperdynamics is enabled, a few extra parameters must be specified in params.in, which holds all bias related parameters for a given simulation run. Hyperdynamics requires (1) a fill limit, (2) definition of initial and product states, (3) the simulation timestep. 

1. Fill limit: This parameter sets the depth of the bias potential. The bias will not fill higher than this limit. Units are kJ/mol.

2. Initial and Product states: Hyperdynamics requires the detection of transition events. Accordingly, the user must define the initial state and its boundaries. The input for the published simulations looks like ```0.9 0.9 1 1``` which specifies that the initial state corresponds to values of CV1 < 0.9 and CV2 < 0.9. The initial state has been exited when either CV1 > 1 **or** CV2 > 1. Thus, the first two entries ```0.9 0.9``` indicate the upper CV1 and CV2 boundaries of the initial state. The second two entries ```1 1``` indicate the lower CV1 and CV2 boundaries for all other regions of the CV space. Units are same as CV units, which are nanometer. A sample hyperdynamics input from our ligand simulations is included in [RUNdirs]/RErun/HYPER-params.in

3. Simulation timestep: This last entry is simply the integration timestep in picosecond units. 

Additionally, we point out that reproducing the published ligand escape time estimates requires two changes, with respect to the free energy simulations, beyond these additional parameter entries. We used a different cylinder position for hyperdynamics and we used minimal backbone restraints of the BRD4 bromodomain. The cylinder inputs are found in [RUNdirs]/RErun/cylpoints-hyperdynamics. The protein restraints are found in [RUNdirs]/RErun/shortbacks.itp. 

**If reproducing our hyperdynamics, remember to use the correct cylpoints file and remember to point the [Rundirs]/RErun/topol.top to the shortbacks.itp file rather than back.itp.**

As explained in the publication (submitted, link to follow), our hyperdynamics is inteded to stop simulation after the initial state is exited so that a new trajectory can begin running in that initial state. Thus, hyperdynamics will terminate when the initial state is exited. A file called *fort.87* will be produced which lists ```time-of-exit trajectory-boost last-instantaneous-boost``` where the first two are most important for comupting mean escapte time and mean boost. The last entry ```last-instantaneous-boost``` can serve as a guide for judging whether or not the fill-depth of the bias is set too high or whether the initial state is inadequately defined.

In [Rundirs]/ALANINE we provide the PBS script that was used to collect 200 reactions in alanine dipeptide simulations.

# Requirements
1. ***Simulation cell*** Currently only cubic, tetragonal and orthorhombic systems are supported (angles = 90 degrees). At this time we do not plan to implement irregular systems. 

2. ***RMSD CVs*** Right now, RMSD is supported so a number of applications should be possible. We will release a distance based CV set soon.

3. ***Initial state cannot be wrapped*** The initial coordiates of the atoms in the CVs cannot be wrapped through the periodic boundaries. We avoid needing to communicate the system topology by satisfying this requirement.

[GROMACS 5.0.5]: https://github.com/gromacs/gromacs/releases/tag/v5.0.5
[VMD]: http://www.ks.uiuc.edu/Research/vmd/
[HERE]: http://scitation.aip.org/content/aip/journal/jcp/143/23/10.1063/1.4937939
[RUNdirs]: https://github.com/BradleyDickson/fABMACS/tree/master/RUNdirs

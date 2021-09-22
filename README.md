# Stokesian dynamics in porous media

#### Software required
  - intel/2018u3
  - python/3.6
  
#### How to use?
  1. The simulation code was written in reduced units in which particle size $r_{p}$, particle diffusive time scale $\tau_{d}$, and thermal energy $k_{B}T$ were used as fundamental units of length, time, and energy, respectively.

  2. To run simulation, first, create a square nanopost array with confinement parameter $\zeta$ and nanopost radius $r_{np}$ <br>
    - Include the required inputs in the input file "post_particle.txt" <br>
    - Compile code to generate executable file and run executable file <br>
    ` ifort module.f90 post_particle_pos.f90 -o post ` <br>
    ` ./post `
    
  3. Compute grand resistance matrix on grids within a unit cell <br>
    - Compile code to generate executable file and run executable file <br>
    ` ifort module.f90 FF_res_mat.f90 -o FF_res_mat ` <br>
    ` ./FF_res_mat `
    
  4. Compute fitting parameters of inverse of pair-wise mobility matrix <br>
    ` python mat_fit.py $r_{p} $r_{np}$ `
    
  5. Run particle trajectory simulation <br>
    - Include the required inputs in the input file "traj.txt" <br>
    - Use main code "traj_2d.f90" for 2D model system and "traj_3d.f90" for 3D model system <br>
    - Compile code to generate executable file and run executable file with trajectory sequence number <br>
    ` ifort module.f90 traj_*.f90 -o traj_* ` <br>
    ` echo "sequence" | ./traj_* `
    
  6. Analysis of particle trajectories (diffusivity, dispersion coefficient, average particle velocity, tortuosity, and cosine of correlation angle) <br>
    - Include the required inputs in the input file "analysis.txt" according to the trajectory data <br>
    - Compile code to generate executable file and run executable file <br>
    ` ifort analysis_module.f90 analysis_main.f90 -o analysis ` <br>
    ` ./analysis`

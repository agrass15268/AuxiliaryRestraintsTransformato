# AuxiliaryRestraintsTransformato
BSc. Thesis: Relative Binding Free Energy Calculations using Auxiliary Restraints in Transformato

This branch contains the underlying data for statistical analysis and plots. Due to size constraints, no trajectories are stored here (or anywhere, unless the CBC doesnt clean them out).

These files were generally created using the tools in [transformato-helpers](https://github.com/agrass15268/transformato-helpers)

### General simulation parameters:

- Hydrogen Mass Repartitioning (either via CHARMMGUI or parmed) was used for all simulations
- Simulation step time: either 0.001 ps for 1.25 ns runs, or 0.004 ps for 5 and 10 ns runs, with the number of steps scaled accordingly
- HBonds were always constrained
- All runs without explicit note used harmonic potentials
- All runs were simulated at 303.15 K

### File types:

- ma_\*.csv files: Contain the results of the mbar analysis, alongside a variety of restraint parameters and the timescale
- relativedistances/\*.csv files: In lieu of trajectory files, these contain csv files of the relative distances relative to the initial distance experienced by the restraints.
  For replicates without restraints, restraints from an ex3 run were used as baseline comparision. Only a subset of trajectory frames was analysed, with a step of 50.

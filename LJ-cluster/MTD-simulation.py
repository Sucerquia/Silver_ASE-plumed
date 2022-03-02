# We used this code to run the LJ cluster benchmark
# simulation using ASE.

from ase.calculators.lj import LennardJones
from ase.calculators.plumed import Plumed
from ase.io.trajectory import Trajectory
from ase.constraints import FixedPlane
from ase.md.langevin import Langevin
from ase.io import read
from ase import units


atoms = read('input.xyz')
timestep = 0.005

i = open("plumed.dat", "r").read().splitlines()

cons = [FixedPlane(i, [0, 0, 1]) for i in range(7)]
atoms.set_constraint(cons)
calc = Plumed(calc=LennardJones(rc=2.5, r0=3.0),
              input=i,
              timestep=timestep,
              atoms=atoms,
              kT=0.1)

atoms.set_masses([1, 1, 1, 1, 1, 1, 1])
atoms.calc = calc
dyn = Langevin(atoms, timestep,
               temperature_K=0.1/units.kB,
               friction=1,
               fixcm=False)

tra = Trajectory('output.traj', 'w', atoms)
dyn.attach(tra.write, interval=100)

dyn.run(1000000)

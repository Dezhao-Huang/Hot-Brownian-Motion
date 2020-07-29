#!/bin/bash
list_lam=" 1000 900 800 700 600 500 400 300 200 90"
#list_lam=" 1000"

for i in $list_lam
do
lam=$i
dir=/scratch365/dhuang2/bubble_2020/1_heat_loop/fix_bubble/${lam}
mkdir -p $dir
cd $dir

cat > heat.in << endofdata

shell rm -rf dump
shell mkdir dump


variable TParticle equal ${lam}

# lammps bubble
units metal        
dimension 3        
boundary p p p     
atom_style atomic  

restart 1000 ./dump/rst.1 ./dump/rst.2

read_restart /scratch365/dhuang2/bubble_2020/0create/large_150/0equil/dump/rst.nve
# # potentials
pair_style hybrid lj/cut 13.5 morse 12
pair_coeff 1 1 lj/cut 	0.0104 3.40
pair_coeff 1 2 lj/cut 	0.02   3.21
pair_coeff 1 3 lj/cut 	0.02   3.21
pair_coeff 2 2 morse  	0.475  1.583  3.0242  12 
pair_coeff 2 3 lj/cut 	0.475  3.0242
pair_coeff 3 3 morse  	0.475  1.583  3.0242  12 

neighbor 2.0 bin
neigh_modify every 10 delay 0 check yes

timestep 0.005
compute ltemp liquid   temp
compute ptemp particle temp
compute wtemp wall      temp
compute nocomtemp particle temp/com

# operations
fix 0  particle spring/self 2.925

fix	1  all nve
fix 2  all press/berendsen  x   1.0 1.0  50.0  y   1.0 1.0  50.0  z   1.0 1.0  50.0

fix 3  particle temp/rescale 1000 \${TParticle} \${TParticle} 1 1
fix_modify   3    temp   nocomtemp

fix 4  wall     langevin 90 90  0.2 699483 tally yes

# output number mass density

variable xcmx equal xcm(particle,x)
variable xcmy equal xcm(particle,y)
variable xcmz equal xcm(particle,z)

compute cc6 liquid chunk/atom bin/sphere \${xcmx} \${xcmy} \${xcmz}   0.0 100.0  40 
fix 611   liquid ave/chunk 1000 10 10000 cc6 temp density/number density/mass file liquid_sphere.profile

# output thermodynamics information
thermo           1000
thermo_style custom step c_ltemp c_ptemp c_wtemp etotal ke pe  press  pxx pyy pzz   lx ly lz 
dump             1 all custom 10000  ./dump/heat.lammpstrj id type x y z vx vy vz fx fy fz 

run				 1000000

undump           1
write_restart    ./dump/rst.fixbubble
write_data       ./dump/data.fixbubble


endofdata


cat > submit.sh << endofdata
#!/bin/csh
#$ -M  notredamecrc@gmail.com    # Email address for job notification
#$ -m  abe        # Send mail when job begins, ends and aborts
#$ -pe smp 64      # Specify parallel environment and legal core size
#$ -q  long@@tengfeiluo
#$ -N  fnb${lam} # Specify job name

module load intel/18.0 ompi/3.0.0/intel/18.0
module load lammps    # Required modules
mpiexec -n \$NSLOTS lmp_mpi < heat.in 
endofdata

qsub submit.sh
cd ..

done

# qstat -u dhuang2
# nodesInUse.sh @tengfeiluo
# echo 'nodesInUse.sh @tengfeiluo'

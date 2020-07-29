#!/bin/bash

temps=" 1000 900 800 700 600 500 400 300 200 90 "
# temps=" 1000 "

for i in $temps
do

temp=$i
mkdir -p  T${temp}
list_lam=" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
 49 50 "
# list_lam=" 1 2 3 "

for j in $list_lam
do

rstnum=$j
dir=/scratch365/dhuang2/bubble_2020/2_diffusion/T${temp}/rst${rstnum}
mkdir -p $dir
cd $dir

cat > heat.in << endofdata

shell rm -rf dump
shell mkdir dump

variable TParticle equal ${temp}


units metal        
dimension 3        
boundary p p p     
atom_style atomic  

restart 1000 rst.1 rst.2

read_restart /scratch365/dhuang2/bubble_2020/1_heat_loop/write_rst/${temp}/dump/rst.${rstnum}
shell echo /scratch365/dhuang2/bubble_2020/1_heat_loop/write_rst/${temp}/dump/rst.${rstnum}

reset_timestep 0

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

fix	1  all nve

fix 3  particle temp/rescale 1 \${TParticle} \${TParticle} 1 1
fix_modify   3    temp   nocomtemp

fix 4  wall     langevin 90 90  0.2 699483 tally yes

# computations

variable xcmx equal xcm(particle,x)
variable xcmy equal xcm(particle,y)
variable xcmz equal xcm(particle,z)
variable vcmx equal vcm(particle,x)
variable vcmy equal vcm(particle,y)
variable vcmz equal vcm(particle,z)
compute myforce particle group/group liquid

compute centermass  particle com
fix centermass particle ave/time 10 1 10 c_centermass[*] file ./dump/com.out mode vector

compute mymsd particle msd
fix mymsd particle ave/time 10 1 10 c_mymsd[*] file ./dump/msd.out mode vector

compute mymsd_L liquid msd
fix mymsd_L liquid ave/time 10 1 10 c_mymsd_L[*] file ./dump/msd_Liquid.out mode vector


variable    p equal 2000     # correlation length
variable    s equal 1       # sample interval
variable    d equal \$p*\$s   # dump interval

fix          SS all ave/correlate \$s \$p \$d &
             v_vcmx v_vcmy v_vcmz c_myforce[1] c_myforce[2] c_myforce[3]  & 
             type auto file ./dump/S0St.dat ave running

variable     v11 equal trap(f_SS[3])
variable     v22 equal trap(f_SS[4])
variable     v33 equal trap(f_SS[5])
variable     v44 equal trap(f_SS[6])
variable     v55 equal trap(f_SS[7])
variable     v66 equal trap(f_SS[8])

variable     vcor equal (v_v11+v_v22+v_v33)/3.0
variable     fcor equal (v_v44+v_v55+v_v66)/3.0


thermo       \$d
thermo_style custom step c_ltemp c_ptemp c_wtemp etotal ke pe  press  pxx pyy pzz lx ly lz v_vcor v_fcor

run				 100000

write_restart    ./dump/free.rst
write_data       ./dump/free.data

endofdata


cat > submit.sh << endofdata
#!/bin/csh
#$ -M  notredamecrc@gmail.com    # Email address for job notification
#$ -m  abe        # Send mail when job begins, ends and aborts
#$ -pe smp 32   # Specify parallel environment and legal core size
#$ -q  long@@tengfeiluo
#$ -N  D${temp}_${rstnum} # Specify job name

module load intel/18.0 ompi/3.0.0/intel/18.0
module load lammps    # Required modules
mpiexec -n \$NSLOTS lmp_mpi < heat.in 
endofdata

qsub submit.sh
cd ..

done
cd /scratch365/dhuang2/bubble_2020/2_diffusion
done
# qstat -u dhuang2
# nodesInUse.sh @tengfeiluo
# echo 'nodesInUse.sh @tengfeiluo'

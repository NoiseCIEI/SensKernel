#!/bin/csh
################################################################################
# C-shell script: KERNELS
################################################################################
# Purpose: to get phase and group velocity sensitivity kernels for the 1-D layered model
################################################################################
# Usage:
#
# % KERNEL model outname wavetype 0-mode End-mode Start-Period End-Period Period-Increment 
################################################################################
set nonomatch
if ( $#argv != 8) then
echo "USAGE: grv_sens_kernel.scr profile outname wavetype(R or L) end_mode start_period end_period incr_period incr_depth"
exit 1 
endif
#-----------------------------
#-----------------------------
#-----------------------------
#-----------------------------
#-----------------------------
#--------Script---------------
set pert=1.0
set out=${2}
../bin/SURF_PERTURB $1 $out $3 0 $4 $5 $6 $7 -s $8 -a -f -p $pert
set out_1=${2}.phv
echo "out_1=$out_1"
../bin/PHV_SENS_KERNEL $1 $out  $3  $out_1 
set out=${2}-
set pert=0.99
../bin/SURF_PERTURB $1 $out $3 0 $4 $5 $6 $7 -s $8 -a -f  -p $pert
set out_1=${2}-.phv
echo "out_1=$out_1"
echo "pert"
../bin/PHV_SENS_KERNEL $1 $out  $3  $out_1
set out=${2}+
set pert=1.01
../bin/SURF_PERTURB $1 $out $3 0 $4   $5 $6 $7 -s $8 -a -f -p $pert
set out_1=${2}+.phv
echo "out_1=$out_1"
../bin/PHV_SENS_KERNEL $1 $out  $3  $out_1
set out_1=${2}.grv
echo "out_1=$out_1"
../bin/GRV_SENS_KERNEL $2 $4  $3  
\rm *-*
\rm *+*
exit
# Example of test runs by script RUN in working directory.
# RUN include two commands:
../bin/KERNELS.csh eus_model test R 1 10 100 5 2
../bin/KERNELS.csh eus_model test L 1 10 100 5 2
#
# eus_model is a standard control profile
#----------------Arguments------
#$1   Name of the file with radial profile of the  elastic and nonelastic parameters (like prem model)
#$2   Name of the output files
#$3   Wavetype (R or L for Rayleigh or Love Waves)
#$4   End-mode (0,1,2,..; f.g., 3). Calculations always start from 0 (fundamental) mode
#$5   Start-period, T0 sec (f.g., 10)
#$6   End-period, Tn sec   (f.g., 150)
#$7   Increment in period, dt sec  (f.g., 10). Number of periods nper=(Tn-T0)/dT+1
#$8   Increment in depth, km (for the fundamental mode 1 or 2 or 3 km are OK; for overtones 10 km is OK)
#----------------Input files-------
#     File with radial profile of the  elastic and nonelastic parameters (like prem model)
#----------------Output files-------
#     File with the layered structure ($2)
#          There are five columns: thickness (km), a (Vp velocity, km/s), b (Vs velocity, km/s),rho (density,g/cm**3),
#          Qs  
#     File with all forward calculation info ($2.wavetype), with complicated structure (may ignore it).
#     File with phase velocity curves for all modes divided by two blank lines 
#         ($2.wavetype.phv)
#         There are 2 columns: period (T), phase velocity (c)
#     File with group velocity curves for all modes divided by two blank lines 
#         ($2.wavetype.grv)
#         There are 2 columns: period (T), group velocity (u)
#     File with apparent Q-factor curves for all modes divided by two blank lines 
#         ($2.wavetype.att)
#         There are 2 columns: period, apparent Q
#     Files with normalized partial derivatives of phase velocity for
#         all requested modes (0,1,..  #     last_mode) and all nper requested periods: 
#         $2.phv.wavetype_mode(j)_period(i) (i=1,2,...,nper),
#         nper is a number of periods.
#         Each file has 3 ($3=L) or 4 ($3=R) columns: depth, [(dc/c)/(db/b)],[(dc/c)/(drho/rho)] or 
#         depth, [(dc/c)/(db/b),[(dc/c)/(da/a)],[(dc/c)/(drho/rho)].
#         The first raw contains also values of period,phase velocity,group velocity,mode_number
#     Files with normalized partial derivatives  of group velocity for 
#         all requested modes (0,1,..,last_#     mode) and all nper requested periods: 
#         $2.grv.wavetype_mode(j)_period(i) (i=1,2,...,nper), 
#         nper is a  number of periods. 
#         Each file has 3 ($3=L) or 4 ($3=R) columns: depth, [(du/u)/(db/b)],[(du/u)/(drho/rho)] or 
#         depth,[(du/u)/(db/b)],[(du/u)/(da/a)],[(du/u)/(drho/rho)].
#         The first raw contains also values of period,phase velocity,group velocity,mode_number
#     General comments:
#         Partial derivatives have dimension 1/km; to calculate actual absolute perturbation in velocity
#         caused by the relative perturbation of parameters a,b,rho between depths Z1 and Z2
#         it is necessary to calculate the corresponding integral by the depth
#         The integrands are:
#         for perturbation of c due to perturbation of b:   [(dc/c)/(db/b)]*c*(del(b)/b);
#         for perturbation of u due to perturbation of b:   [(du/u)/(db/b)]*u*(del(b)/b);
#         for perturbation of c due to perturbation of a:   [(dc/c)/(da/a)]*c*(del(a)/a);
#         for perturbation of u due to perturbation of a:   [(du/u)/(da/a)]*u*(del(a)/a);
#         for perturbation of c due to perturbation of rho: [(dc/c)/(drho/rho)]*c*(del(rho)/rho);
#         for perturbation of u due to perturbation of rho: [(du/u)/(drho/rho)]*u*(del(rho)/rho).
#         Here:
#          del(b) is actual perturbation of b, (km/s);
#          del(a) is actual perturbation of a, (km/s);
#          del(rho) is actual perturbation of rho, (g/cm**3).

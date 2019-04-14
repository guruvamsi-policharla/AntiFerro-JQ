#!/bin/bash
##### Queue #####
#PBS -q Small
###### Select resources #####
#PBS -N driver_mc_4
#PBS -l nodes=1:ppn=10,walltime=50:00:00
#PBS -l mem=4g
###### Priority #####
#### Output File #####
#PBS -o $PBS_O_WORKDIR/4_2mc.out
#### Error File #####
#PBS -e $PBS_O_WORKDIR/4_2mc.err
##### Change to current working directory ##### cd $PBS_O_WORKDIR
##### Execute Program #####
hostname
julia /home/guru/repos/antiFerro/main.jl
date

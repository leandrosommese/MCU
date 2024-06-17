#!/bin/bash
pdb_basename=$1
ff=$2


## Select the Force Field: From '/usr/local/gromacs/share/gromacs/top'
gmx pdb2gmx -f ${pdb_basename}_clean.pdb -o ${pdb_basename}_processed.gro -water spce -ff ${ff}

# Define box
gmx editconf -f ${pdb_basename}_processed.gro -o ${pdb_basename}_newbox.gro -c -d 1.0 -bt cubic

# Solvation
gmx solvate -cp ${pdb_basename}_newbox.gro -cs spc216.gro -o ${pdb_basename}_solv.gro -p topol.top

# Add ions
## Requires ions.mdp file
gmx grompp -f ions.mdp -c ${pdb_basename}_solv.gro -p topol.top -o ions.tpr

# Doc [https://manual.gromacs.org/current/onlinehelp/gmx-genion.html]
echo ${ion_replace_group} | gmx genion -s ions.tpr -o ${pdb_basename}_solv_ions.gro -p topol.top -np 0 -pname NA -pq 1 -nn 0 -nname CL -nq -1 -conc 0 -neutral

# Energy minimization
## Requires minim.mdp file
gmx grompp -f minim.mdp -c ${pdb_basename}_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
echo ${potential_data} | gmx energy -f em.edr -o potential.xvg

# Equilibration
## Requires nvt.mdp and npt.mdp files
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

echo ${temperature_data} | gmx energy -f nvt.edr -o temperature.xvg # temperature progression

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt 

echo ${pressure_data} | gmx energy -f npt.edr -o pressure.xvg # pressure progression

echo ${density_data} | gmx energy -f npt.edr -o density.xvg # density progression

# Production MD
## Requires md.mdp file
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1

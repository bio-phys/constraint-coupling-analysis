# Constraint coupling analysis

Compute the largest absolute eigenvalue of the A matrix in the [LINCS](https://doi.org/10.1021/ct700200b) algorithm for the [standard](https://doi.org/10.1021/jp207925m), [HMR](https://doi.org/10.1021/acs.jctc.9b00160) and [VIS](https://doi.org/10.1021/acs.jctc.8b00267) CHARMM versions. This repository only contains a short script and example Gromacs files. No installation necessary.

## Requirements:
* python3
* [numpy](https://numpy.org/)
* [MDAnalysis](https://www.mdanalysis.org/)

## How-to:
1) Depending on your [Gromacs](https://www.gromacs.org/) and MDAnalysis versions, you might have to re-create the tpr files:
    ``` console
    gmx grompp -f standard.mdp -c standard.gro -p standard.top -o standard.tpr
    gmx grompp -f hmr.mdp -c hmr.gro -p hmr.top -o hmr.tpr
    gmx grompp -f vis.mdp -c vis.gro -p vis.top -o vis.tpr
    ```

2) MDAnalysis does not have a separate representation for bonds and constraints. You need to specify a file containing the indices of the bond/constraint to be included in the analysis. You can get a list of constraints from
    * the itp file
    * the tpr file (using `gmx dump -s *.tpr`)
    * an interactive / debugging MDAnalysis session

    For ease of use, the list of selected constraints is printed to screen.
    **NOTE:** the list must be indexed from 0 (see for example `constr-all.dat`)

3) Run the analysis:
    ``` console
     python constraint-coupling.py datafiles/hmr.tpr datafiles/hmr.gro constr-all.dat
    ```
    **NOTE:** a negative LINCS order ususally means that the eigenvalue is larger than 1.


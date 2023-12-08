# rcladder-id

## Structured Identification of RC-Ladder-Type Systems
A set of MATLAB scripts to perform structured system identification of SISO RC-ladder type systems of the form:
![RC-ladder schematic](/img/rcladder.PNG)

The code is divided into several parts:

## Main Scripts

| Name | Description | Link to paper |
|------|-------------|---------------|
| `IDMonteCarlo.m` | Main script to run the algorithm for the examples in [1] | Section III |
| `GenerateFigures.m` | Script to visualize data from the paper | Figure 3 |

### Helper Functions

| Name | Description | Link to paper |
|------|-------------|---------------|
| `RCLadderN.m` | Generates an admittance state-space model from a vector of R and C parameters. | Equation 2 |
| `RCLadder2Theta.m` | Performs the conversion between a state-space model already placed in the correct form into a vector of resistance and capacitance parameters (inverse of `RCLadderN.m`). | Equation 2, line (7) from Algorithm 1. |
| `RCLadderDiagScaling.m` | Performs the procedure to restore the unknown diagonal scaling. | Section II-B, line (10) from Algorithm 1. |

### Tridiagonal Matrix Operations

The folder `/Tridiagonal Matrix Operations` containing script to manipulate tridiagonal matrices:

| Name | Description | Link to paper |
|------|-------------|---------------|
| `tridiagcrout.m` | Crout LU factorization | Section II-E, line (20) from Algorithm 1. |
| `tridiaglanczos.m` | Performs the orthogonal Lanczos process | Section II-B, line (5) from Algorithm 1. |
| `tridiagpartlett.m` | Partlett tridiagonalization procedure. | Not used. |
| `tridiagsim.m` | Tridiagonalization via similarity transformations. | Not used. |
| `tridiagsymmetrize.m` | Compute the similarity transformation to symmetrize a tridiagonal matrix. | Not used. |

### Other Methods

The `/Other Methods` folder contains implementations of:

1. `/Calzavara2021`: reference [9] from the paper (partial implementation).
2. `/Hwang1984`: reference [7] from the paper.
3. `/Yu2018`: reference [18] from the paper.

The methods 2. and 3. from above were used in the comparison from the paper.

### img

This folder includes graphics for Figure 3.

### Data

This folder includes data from Figure 3.

## More Information

Please cite if this work is useful to you. These scripts were used in the following publication:

> Hannigan, B. C. & Menon, C. (2023). Fast, analytical method for structured identification of SISO RC-ladder-type systems. _IEEE Transactions on Circuits and Systems II: Express Briefs_, 1. https://doi.org/10.1109/TCSII.2023.3340505

©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon
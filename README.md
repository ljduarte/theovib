# Theovib
[![GitHub issues](https://img.shields.io/github/issues/ljduarte/theovib)](https://github.com/ljduarte/theovib/issues)
[![GitHub forks](https://img.shields.io/github/forks/ljduarte/theovib)](https://github.com/ljduarte/theovib/network)
[![GitHub stars](https://img.shields.io/github/stars/ljduarte/theovib)](https://github.com/ljduarte/theovib/stargazers)
[![GitHub license](https://img.shields.io/github/license/ljduarte/theovib)](https://github.com/ljduarte/theovib/blob/main/LICENSE)

:construction: Version 0.0.1 :construction:

A package for molecular vibrations analysis. This initial version contains functions to solve the vibrational problem starting from the Hessian matrix or from a 3D-Hessian construct from Interacting Quantum Atoms energy decomposition scheme. 
Infrared intensities are obtained from atomic charges and dipoles obtained by AIMAll.   

## Contents
[Installation](#Installation)

[Usage](#Usage)

[Theory](#theory)

[Contributing](#Contributing)

[Documentation](#Documentation)

[License](#License)

## Installation

To install, just use pip:

```bash
$ pip install theovib
```

## Usage

You can perform a vibrational analysis for water using the files in the **"example"** folder.

### Vibrational analysis of the Water molecule

In **./example/h2o** you will find all single point calculations needed to generate the 3D-Hessian matrix.

The Cartesian coordinates of the water molecule are:

| **Atom(Label)** |    **X** |     **Y** |     **Z** |
|-----------------|---------:|----------:|----------:|
| O(1)            | 0.000000 |  0.000000 |  0.004316 |
| H(2)            | 0.000000 | -0.763369 | -0.580667 |
| H(3)            | 0.000000 |  0.763369 | -0.580667 |

The internal coordinates are defined as:

* **bond** between atoms **1** and **2**;
* **bond** between atoms **1** and **3**;
* **angle** defined by atoms **1**, **2** and **3**.

We start writing the input file:

```bash
MOLECULE:
h2o
FOLDER:
h2o
DELTA:
0.05
BOND:
1 2
1 3
---
ANGLE:
2 1 3
---
```

**DELTA** is the displacement that generated the non-equilibrium geometries.

#### 1. Reading the input file:

Import modules:

```python
 from theovib.molecule import *
 from theovib.internal import *
 from theovib.matrices import *
 from theovib.ir import *
 from theovib.input import *
```

Use the **Input** class:

```python
input_data = Input.read_text('input.txt')
```

Use the **Molecule** class and its methods to get and store data:

```python
molecule = Molecule.read_gaussian(input_data.folder + '/EQ.com')
molecule.energy = get_energy_from_wfn(input_data.folder +'/EQ.wfn')
molecule.iqa_energy = get_IQA(input_data.folder +'/EQ_atomicfiles', molecule.atoms)
```

#### 2. Construct the B matrix:

Initialize and construct the **B** matrix:

```python
b_matrix = []
    for coord in input_data.bond:
        b_matrix.append(bond(molecule.positions, coord[0], coord[1]))
    for coord in input_data.angle:
        b_matrix.append(angle(molecule.positions, coord[0], coord[1], coord[2]))
molecule.b_matrix = np.array(b_matrix)    
```

#### 3. Obtain the Hessian and 3D Hessian matrix:

```python   
molecule.hessian, molecule.iqa_hessian, errors = hessian_from_iqa(molecule.atoms, input_data.delta, input_data.folder)
```

The Hessian and 3D Hessian are numerrically generated using IQA contributions from AIMAll outputs

#### 3. Calculate the normal coordinates and infrared intensities:

```python
molecule.normal_coordinates, molecule.freq, molecule.iqa_freq, molecule.iqa_terms = normal_modes(molecule.atoms, molecule.iqa_hessian)
molecule.int, molecule.c_tensors, molecule.ct_tensors, molecule.dp_tensors = intensities(molecule.atoms, molecule.positions, molecule.normal_coordinates, input_data.folder, input_data.delta)
```

#### 4. Convert to internal coordinates to obtain the force constants:

```python
molecule.internal_hessian, molecule.iqa_forces = convert_to_internal(molecule.atoms, molecule.b_matrix, molecule.iqa_hessian)
```

## Theory
#### 1. Force constant in internal coordinates:

In order to calculate, and decompose the force constants into its IQA components, one needs first to convert the Cartesian Hessian into the the Wilson's $\mathbf{F}$ matrix, that contains the force constants, and their interactions, expressed in internal coordinates. To do this, it is necessary to define the $\mathbf{B}$ matrix, that converts the $N\times 3N$ Cartesian coordinates matrix, $\mathbf{X}$ in the internal coordinates matrix $\mathbf{R}$:

$$ \mathbf{B}\mathbf{X} = \mathbf{R} $$

The process of setting-up the $\mathbf{B}$ matrix is tedious and can be found in the literature. 
We start by calculating the pseudo-inverse of $\mathbf{B}$ by:

$$ \mathbf{B}^{-1} =\mathbf{M}^{2} \mathbf{B}^{\dagger}\mathbf{G}^{-1} $$

$\mathbf{G}$ contains the inverse of the kinetic energy terms and its inverse, $\mathbf{G}^{-1}$, is given by:

$$\mathbf{G}^{-1} =\mathbf{D} \mathbf{\Phi}^{-1}\mathbf{D}^{\dagger}$$

where $\mathbf{D}$ and $\mathbf{\Phi}$ are, respectively, the eigenvectors and the diagonal eigenvalues matrices of $\mathbf{G}$. The force constants in internal coordinates are then obtained as follows:

$$\mathbf{F} ={\mathbf{B}^{\dagger}}^{-1} \mathbf{H}\mathbf{B}^{-1}$$   

The decomposition of the force constants into the IQA contributions is done using:

$$\left[\sum_{k=1}^{N^2} \mathbf{F^{IQA}_k}\right] ={\mathbf{B}^{\dagger}}^{-1} \left[\sum_{k=1}^{N^2} \mathbf{H^{IQA}_k}\right] \mathbf{B}^{-1}$$

Each $F^{IQA}_k$ is a matrix containing the contribution of the $k^{th}$ IQA term of $\mathbf{F}$. 

#### 2. Infrared intensities:

The infrared intensity, $A_k$, of normal mode $Q_k$ is defined as:

$$A_k = \frac{N_A \pi}{3c^2}\left(\frac{\mathrm{d}\vec{p}}{\mathrm{d}Q_k}\right)^2 = \frac{N_A \pi}{3c^2} \sum_{i=1}^3 \left( \frac{\partial\vec{p}}{\partial \sigma_i} \cdot \frac{\partial\sigma_i}{\partial Q_k}\right)^2$$

where $N_A$ is the Avogadro's constant, $c$ is the speed of light and $\vec{p}$ is the molecular dipole moment. The derivative $\frac{\partial\sigma_i}{\partial Q_k}$ is a element of the $k^{th}$ column of $\mathbf{L}$.

$$\mathbf{L} = \mathbf{M}\mathbf{A}$$

Let $\mathbf{T}$  be a block-diagonal matrix of dimension $3N \times 3N$ where each $3 \times 3$ block is an atomic polar tensor. The dipole moment derivative of normal mode $Q_k$ is obtained by solving:

$$\frac{\mathrm{d}p}{\mathrm{d}Q_k} = \mathbf{U} \cdot \mathbf{T} \cdot \mathbf{L_k}$$
$\mathbf{L_k}$ is the $k^{th}$ column of $\mathbf{L}$. $\mathbf{U}$ is a $1 \times N$  line vector whose elements equals 1. 

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## Documentation
https://ljduarte.github.io/theovib/

## License

`theovib` was created by L. J. Duarte. It is licensed under the terms of the MIT license.

## Credits

`theovib` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).

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

[Contributing](#Contributing)

[License](#License)


## Installation

To install, just use pip:

```bash
$ pip install theovib
```

## Usage
You can perform a vibrational analysis for water using the files in the "example" folder.

### Vibrational analysis of the Water molecule
In ./example/h2o you will find all single point calculations needed to generate the 3D-Hessian matrix.

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

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`theovib` was created by L. J. Duarte. It is licensed under the terms of the MIT license.

## Credits

`theovib` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).

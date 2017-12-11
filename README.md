ChemKinLib [![Build Status](https://travis-ci.org/cs207-group11/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/cs207-group11/cs207-FinalProject.svg?branch=master) [![Coverage Status](https://coveralls.io/repos/github/cs207-group11/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/cs207-group11/cs207-FinalProject?branch=master)
===================

### Group #11: Karan R. Motwani, Shiyu Huang, Hannah Sim, and Haixing Yin.

ChemKinLib (Chemical Kinetics Library) consists of classes and functions to allow the user to interpret the kinetics of a system of chemical reaction by computing the associated reaction rates. The library currently supports:

- Elementary reactions: irreversible and reversible

And of these reactions, the following reaction rate coefficients are supported:
- Constant
- Arrhenius
- Modified Arrhenius

This Python library is designed to exhibit flexibility, extensibility, and ease of use. The initial input provided by the user is an `.xml` file while other individual parameters (temperature and concentrations) are entered during runtime.

## Installation (to be updated)
1) Clone the repository into you local system.

2) Open the 'parser.py' file and changed the variable 'path' to the location of your XML file (to be improved using `argparse`).

3) Run chemkin.py via terminal.

4) Parameters specific to each reaction read from your XML file will be requested from you at run-time. (Ex. Concentrations, Temperature)

5) Reaction rates for each individual reaction are then displayed on the screen.

## Documentation
The documentation for the library can be found [here](https://github.com/cs207-group11/cs207-FinalProject).

## Dependencies
- numpy
- scipy
- pytest

## License
See [LICENSE](https://github.com/cs207-group11/cs207-FinalProject/blob/master/LICENSE) file distributed with ChemKinLib for more information.

## Contributing
Contributions are welcome (post project evaluation). If you wish to contribute, please take a few moment to review the [branching model](http://nvie.com/posts/a-successful-git-branching-model/) this repository utilizes.

## Support
If you have questions or need help with using or contributing to the repository, feel free to ask questions through:
- https://github.com/cs207-group11/cs207-FinalProject/issues

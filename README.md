`ChemKinLib` [![Build Status](https://travis-ci.org/cs207-group11/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/cs207-group11/cs207-FinalProject.svg?branch=master) [![Coverage Status](https://coveralls.io/repos/github/cs207-group11/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/cs207-group11/cs207-FinalProject?branch=master)
===================

### Group #11 Members: Karan R. Motwani, Shiyu Huang, Hannah Sim and Haixing Yin.

ChemKinLib stands for Chemical Kinetics Library.

ChemKinLib consists of algorithms, functions, and documentation to allow the user to interpret the kinetics of a chemical reaction by returning the reaction rate for a system of chemical reaction(s). The library currently supports:

- Elementary reactions
- Irreversible reactions
- Reactions with constant reaction rate coefficients
- Reactions with Arrhenius reaction rate coefficients
- Reactions with modified Arrhenius reaction rate coefficients

This Python library is designed to have flexibility, extensibility, and ease of use. Our software design follows an object-oriented approach. The initial input provided by the user is an XML file while the individual parameters for a specific reaction are entered during runtime.

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
- pytest

[Note : Anaconda's standard Python platform should suffice the dependency requirements.]

## License
See [LICENSE](https://github.com/cs207-group11/cs207-FinalProject/blob/master/LICENSE) file distributed with ChemKinLib for more information.

## Contributing
Contributions are welcome (post project evaluation!). If you wish to contribute, please take a few moment to review the [branching model](http://nvie.com/posts/a-successful-git-branching-model/) this repository utilizes.

## Support
If you have questions or need help with using or contributing to the repository, feel free to ask questions through:
- https://github.com/cs207-group11/cs207-FinalProject/issues

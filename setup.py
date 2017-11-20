
import os
from setuptools import setup, find_packages

# This reads the __version__ variable from chemkinlib/_version.py
exec(open('src/chemkinlib/_version.py').read())

# README file as long_description:
long_description = open('README.md').read()

# Read in requirements.txt
requirements = open('requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

setup(
    name='chemkinlib11',
    version=__version__,
    description='Chemical kinetics code',
    long_description=long_description,
    install_requires=requirements,
    url='https://github.com/cs207-group11/cs207-FinalProject',
    author='ChemKinLib Developers',
    author_email='hsim13372@gmail.com',
    license='MIT',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    zip_safe=False,
    include_package_data=True,
    package_data={
        '': [os.path.join('src', 'chemkinlib', 'data', '*.xml'),
             os.path.join('src', 'chemkinlib', 'data', '*.sqlite')]
    }
    )
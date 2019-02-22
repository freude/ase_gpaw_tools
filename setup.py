from setuptools import setup, find_packages


with open("README.md", "r") as fh:
        long_description = fh.read()

setup(name='ase-gpaw-tools',
      version='0.1',
      description='Library of scripts facilitating ASE-GPAW computations',
      author='Mike Klymenko',
      author_email='mike.klymenko@rmit.edu.au',
      license='MIT',
      packages=find_packages(),
      zip_safe=False
      )


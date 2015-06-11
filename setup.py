# setuptools installation of proto_md
# Copyright (c) 2011-2015 Adrew Abi-Mansour : anabiman at indiana dot edu
# Released under the GNU Public License 2 (or higher, your choice)
# setup.py originally boosted from Oliver Beckstein's Gromacs Wrapper 
#
from __future__ import with_statement

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

with open("README") as readme:
    long_description = readme.read()


version = "0.0.1"

setup(name="DMS",
      version=version,
      description="A library for doing coarse grained brownian motion dynamics.",
      long_description=long_description,
      author="Andy Somogyi :: Andrew Abi-Mansour",
      author_email="anabiman at indiana dot edu",
      license="GPLv2",
      url="://github.com/CTCNano/proto_md",
      keywords="science Gromacs analysis 'molecular dynamics'",
      classifiers=['Development Status :: 0.0.1 - way pre alpha',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Operating System :: POSIX',
                   'Operating System :: Linux',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ],
      packages=find_packages(exclude=['tests','scripts','extras','doc/examples']),
      scripts = [],
      package_data={'proto_md': ['templates/*.sge', 'templates/*.pbs',  # template files
                                'templates/*.ll', 'templates/*.sh',
                                'templates/*.mdp', 'templates/*.cfg',
                                'external/GridMAT-MD_v1.0.2/GridMAT-MD.pl']   # external bundled scripts
                    },
      install_requires = ['numpy>=1.0',
                          'scipy',        # numkit needs it
                          "MDAnalysis", 
                          "GromacsWrapper"
                          ],              # basic package (w/o analysis)
      extras_require = {
                'numkit': ['scipy'],
                },
      zip_safe = True,
)

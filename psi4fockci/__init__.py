#
# @BEGIN LICENSE
#
# psi4fockci by Shannon E. Houck, a plugin to perform 
# generalized Fock-space CI calculations using the 
# DETCI module.
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""
A program for running RAS-SF-IP/EA using Psi4.


This program runs RAS-nSF-IP/EA using Psi4's DETCI module. The 
program is primarily run through the sf_cas function call.
"""
__version__ = '0.1'
__author__  = 'Shannon E. Houck'

# Load Python modules
from .spinflip import *

# Load C++ plugin
#import os
#import psi4
#plugdir = os.path.split(os.path.abspath(__file__))[0]
#sofile = plugdir + '/' + os.path.split(plugdir)[1] + '.so'
#psi4.core.plugin_load(sofile)

psi4.driver.procedures['energy']['psi4fockci'] = run_psi4fockci


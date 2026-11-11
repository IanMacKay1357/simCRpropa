from crpropa import *
import crpropa
import logging
import yaml
import numpy as np
import argparse
from os import path, environ
from copy import deepcopy
from glob import glob
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
import simCRpropa
import socket
import h5py
import os
from numpy.random import RandomState

def init_custom_rectangular_prism_bfield(
    void_loc,
    vgrid,
    obsSize,
    b_void=0 * gauss,
    b_ext=1e-13 * gauss,
    seed=0,
):
    """Fill grid with magnetic field values where user specifies void location as a rectangular prism shell.
    PUT ALL UNITS IN SI.

    Parameters
    ----------
    void_loc : list of 2 floats
        First number is fraction of LOS where void starts. Second number is where it ends. Creates a rectangular shell
        and sets B-field to b_void within. Sets b_ext != 0 everywhere else.
    b_void : float
        Magnetic field in voids, in Gauss.
    b_ext : float
        Magnetic field outside/exterior to voids, in Gauss.
    seed : int
        Seed for random direction of magnetic field at each grid point.

    Returns
    -------
    """
    print(f"Using random seed={seed}")
    prng = RandomState(seed)
    gridArray = vgrid.getGrid()
    nx = vgrid.getNx()
    ny = vgrid.getNy()
    nz = vgrid.getNz()
    msg = f"vgrid: nx={nx}, ny={ny}, nz={nz}"
    print(msg)

    # Fill grid
    for xi in range(0, nx):
        for yi in range(0, ny):
            for zi in range(0, nz):
                vect3d = vgrid.get(xi, yi, zi)
                x = prng.uniform(-1, 1)
                y = prng.uniform(-1, 1)
                z = prng.uniform(-1, 1)
                d = np.sqrt(x * x + y * y + z * z)

                # if B ≠ 0 at the source, I need to make B = 0 beyond the source to avoid interpolation errors described below.
                if void_loc[0] < 1 <= void_loc[1]:
                    # Check if grid x coordinate is in the user specified void location along the LoS from observer to source (x is along LoS: src->obs)
                    # If B ≠ 0 near the source, pad edge of grid with B=0 so that the center of the grid (the source) does not have B=/=0 from interpolating the edge.
                    # In other words, add B = 0 outside source to avoid edge effects in which the center of the grid has B =/=0 (interpolated from case when B=/=0 at the edge of the grid?)
                    if int(obsSize)*void_loc[0] <= xi <= int(obsSize)*void_loc[1] or xi > obsSize:
                        vect3d.x = b_void * x / d
                        vect3d.y = b_void * y / d
                        vect3d.z = b_void * z / d
                    else:
                        vect3d.x = b_ext * x / d
                        vect3d.y = b_ext * y / d
                        vect3d.z = b_ext * z / d

                else:
                    # B is nonzero near source
                    # It should not matter that B=/=0 beyond the observer (near the grid edge) in this case, because if interpolated
                    # to the center, the desired behavior is B=/=0 there anyway
                    if int(obsSize)*void_loc[0] <= xi <= int(obsSize)*void_loc[1]: #or xi > obsSize:
                        vect3d.x = b_void * x / d
                        vect3d.y = b_void * y / d
                        vect3d.z = b_void * z / d
                    else:
                        vect3d.x = b_ext * x / d
                        vect3d.y = b_ext * y / d
                        vect3d.z = b_ext * z / d
    msg = 'Grid done'
    print(msg)

    return None

def test(config_path):

    # Number of cells to pad grid to avoid edge effects. Default should be 2.
    gridPad = 2
    gridOrigin = Vector3d(0, 0, 0)

    with open(config_path, "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    gridSpacingMpcImplicit = config['Bfield']['geo']['gridSpacing']
    gridSpacingMeters = gridSpacingMpcImplicit * Mpc
    obsSize = int(np.ceil(redshift2ComovingDistance(config['Source']['z']) / Mpc / gridSpacingMpcImplicit ))
    gridSize = obsSize + gridPad
    if gridSize % 2 != 0:
        print(f"The grid size {gridSize} does not have an even length to divide in half. Increasing the grid size by 1.")
        gridSize += 1
    vgrid = Grid3f(gridOrigin,
            gridSize,
            gridSpacingMeters)

    init_custom_rectangular_prism_bfield(void_loc=config['Bfield']['geo']['void_loc'],vgrid=vgrid,obsSize=obsSize,b_void=config['Bfield']['b_void'],b_ext=config['Bfield']['geo']['B'],seed=config['Bfield']['geo']['seed'])

    return vgrid

def plot_vgrid(vgrid):
    
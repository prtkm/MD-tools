
# TODO

# Slope Optimizer
# Read QE output
# Plot thermostat
# Plot conductivities
# Plot ln(D) vs 1/T
# Get slope

import numpy as np
import matplotlib.pyplot as plt
from pycse import regress
try:
    import pymatgen as mg
    from pymatgen.io.zeoio import *
    from pymatgen.io.aseio import AseAtomsAdaptor as aseio
    from pymatgen.transformations.standard_transformations import OxidationStateRemovalTransformation as remox
except:
    pass
from ase.io import read
from subprocess import Popen, PIPE


def cif2cssr(cif, remove = ['Li']):
    '''
    Converts cif files to Zeo++ CSSR files, removes oxidation states and species specified in remove. Must have pymatgen installed. Structure must be ordered
    '''
    filename  = cif.split('.')[-2]
    s = mg.read_structure(cif)
    rem = remox()
    s = rem.apply_transformation(s)
    if remove != None:
        s.remove_species(remove)

    try:
        cssr = ZeoCssr(s)
        cssr.write_file('{0}.cssr'.format(filename))
    except:
        cssr = None
    return cssr

def find_channels(file, probe_radius = 0.5, accuracy = 'normal', rad_file = None):

    '''
    Use zeo++ to find conducting channels in given structure- cif or cssr. Must specify zeo executable as $ZEO in .bashrc for this to work.
    '''
    filename, ext = file.split('.')

    if ext == 'cif':

        cssr = cif2cssr(file)

        if cssr == None:

            return None, None

        file = '{0}.cssr'.format(filename)


    if accuracy !='high':
        if rad_file == None:
            cmd = '$ZEO -chan {0} {1}'.format(probe_radius, file)
        else:
            cmd = '$ZEO -r {2} -chan {0} {1}'.format(probe_radius, file, rad_file)
    else:
        if rad_file == None:
            cmd = '$ZEO -ha -chan {0} {1}'.format(probe_radius, file)
        else:
            cmd = '$ZEO -ha -r {2} -chan {0} {1}'.format(probe_radius, file, rad_file)

    p = Popen(cmd, shell =True, stdout =PIPE, stderr= PIPE)
    out, err = p.communicate()
    if err !='':
        raise Exception('\n\n{0}'.format(err))

    f = open('{0}.chan'.format(filename), 'r')
    lines = f.readlines()

    channels =  lines[0].split()[1]
    dimensionality = lines[0].split()[-1]
    if int(channels) == 0:
        dimensionality = 0

    return int(channels), int(dimensionality)

def read_msd(filepath, skiprows = 1):
    '''
    This function reads a typical msd output file and returns t, msd
    Args : filepath - name of the msd file with two columns t and msd
           skiprows - number of rows to skip on top of the msd file
    Returns: t, msd
    '''

    t, msd = np.loadtxt(filepath, skiprows = skiprows, unpack = True)
    return t, msd



def plot_msd(save = False, show = False, t = None, msd = None, t_units = 'ps', msd_units = '$\AA^{2}$', file = None, skiprows = 1, slope = True, T = None, f = 0.5):
    # TODO make this simpler!!!
    '''
    Plots a mean square displacement trajectory
    Args: save = filepath to save
          show = Turns off/on plt.show()
          Note: Both savefig and show can be done in the main function also
          file = file to read from if t and msd are not specified
          slope = Turns slope on/off
          T = Temperature for legend
          f = fraction of data to discard for equilibration
    '''

    if (t == None or msd == None):
        try:
            t, msd = read_msd(file, skiprows = skiprows)
        except:
            print "Enter either t,msd or file to read from"

    if T != None:
        plt.plot(t, msd, label = '{0} K'.format(T))
    else:
        plt.plot(t,msd)
    if slope == True:
        # Cutting out the equilibration data and using the rest to fit slope
        t_cut = t[len(t)*f:]
        msd_cut= msd[len(t)*f:]

        # Fitting using linear regression and 95 percent confidence intervals - p = [slope, intercept]

        t_stack = np.column_stack([t_cut**1, t_cut**0])
        p, pint, se = regress(t_stack, msd_cut, 0.05)
        msd_fit = np.dot(t_stack,p)

        plt.plot(t_cut, msd_fit)


    plt.xlabel('Time ({0})'.format(t_units))
    plt.ylabel('MSD ({0})'.format(msd_units))

    if T !=None:
        plt.legend()

    if save != False:
        plt.savefig(save)

    if show == True:
        plt.show()

    if slope == True:
        return p, pint, se


def get_D(slope, interval = None):

    D = slope/6. * 1e-4 # in cm^2/s

    if interval == None:
        return D

    else:
        Dint =  np.array(int)/6. * 1e-4
        return D, Dint

def get_conductivity(atoms, T, D = None, slope = None, interval = None,  species = 'all'):

    if D == None:
        if interval == None:
            D = get_D(slope)
        else:
            D, Dint = get_D(slope, interval)

    # Calculating Sigma

    q = 1.60e-19 # Coulombs

    kb = 1.3806488e-23 # Boltzmann Constant in SI units

    if species !='all':
        N = len([atom for atom in atoms if atom.symbol == species])
    else:
        N = len(atoms)

    V = atoms.get_volume() * 1e-24 # cell volume in cm^{3}

    prefactor = N/V * (q**2) /kb / T

    sigma = prefactor * D

    if interval == None:
        return sigma

    else:
        sigma_int = prefactor * Dint
        return sigma, sigma_int

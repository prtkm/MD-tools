
# TODO

# Slope Optimizer
# Read QE output
# Plot thermostat

import numpy as np
import matplotlib.pyplot as plt
from pycse import regress

import pymatgen as mg
from pymatgen.io.zeoio import *
from pymatgen.io.aseio import AseAtomsAdaptor as aseio

from ase.io import read
from subprocess import Popen, PIPE
from pymatgen.transformations.standard_transformations import AutoOxiStateDecorationTransformation as oxi
import os
import shutil

def get_cif_files(dir):

    '''
    Returns a list of all the cif files in a directory
    '''
    files = [file for file in os.listdir(dir) if file.endswith('cif')]
    return files

def remove_duplicates(files, separator):

    '''
    Takes a list of cif files (with extension) and filters out the duplicates.
    File should be of the form <ID><separator><formula>. Usually required while analysing
    ICSD structures where different ids might have the same structure.
    
    Reutrns dict of structure names mapped to a single file.
    This just assigns the last found structure to the file name for now.

    Warning: Matching is only done by formula name. This will also remove all polymorphs of the structure with the same formula, but supercells with unreduced formula wont be removed. Handle with care! 
    '''

    names = [file.strip('.cif').split(separator)[-1] for file in files]
    d = dict(zip(names, files))
    return d


def make_fileset(filenames, dir1, dir2):

    '''
    The point of this function is to make a subset of files (usually cif files for geometry analysis by zeo++).

    Selects files specified in filenames from dir1 and dumps it into dir2. dir1 and dir2 are relative
    to the current directory, unless an absolute path is used.
    '''

    import shutil

    
    if not os.path.isdir(dir2):
        os.makedirs(dir2)

    if not dir1.endswith('/'):
        dir1+='/'
    for file in filenames:
        filepath = dir1+file
        shutil.copy(filepath, dir2)


def prepare_for_md(dir):

    '''
    Filters out cif files with disorder.
    Creates new directories:
    ./md-ready
    ./disordered
    '''
    if not dir.endswith('/'):
        dir+='/'

    files = get_cif_files(dir)
    ordered, disordered = [], []

    for file in files:

        s = mg.read_structure(dir+file)
        if s.is_ordered:
            ordered.append(file)
        else:
            disordered.append(file)

    for subd in ['md-ready', 'disordered']:

        # This creates subdirectories if not already present
        if not os.path.isdir('{0}/{1}'.format(dir, subd)):
            os.makedirs('{0}/{1}'.format(dir, subd))

    for file in ordered:
        shutil.copy(dir + file, dir + 'md-ready/')

    for file in disordered:
        shutil.copy(dir + file, dir + 'disordered/')
        
    return
    
def prepare_for_zeo(dir, remove_duplicates = False, separator = None):

    '''
    Takes a directory of cif files and filters zeo++ usable files. Requires pymatgen.
    Creates new directories:
    ./ready - oxidation state decorated files for zeo++ to use.  
    ./no-rad - structures with species having no corresponding ionic radii in pymatgen database
    ./fails - structures which cannot be assigned oxidation states

    If remove_duplicates is set to true, it runs the files through remove_duplicates. Files pulled out of the database are of the format <id><speparator><formula>. The separator is usually '_' or '-'.
    '''

    # TODO maybe let the user specify the location of the new directories
    
    files = [file for file in os.listdir(dir) if file.endswith('cif')]
    
 
    d = {} # Dictionary of the form d[<filename>]['struct'], d[<filename>]['mass'], d[<filename>]['radius']
           
    fails, no_rad = [], []

    for file in files:
        d[file] = {}
        # reading to pymatgen structure object
        s = mg.read_structure('{1}/{0}'.format(file, dir))
        
        # AutoOxidationStateDecorationObject
        ox = oxi()

        try:
            # Oxidation state decorated structure
            s_ox = ox.apply_transformation(s)

            # Saving structure to dictionary
            d[file]['struct']  = s_ox

            # List of unique elements in the structure
            species = set(s_ox.species)

            radii = dict((str(sp), float(sp.ionic_radius)) for sp in species)
            masses = dict((str(sp), float(sp.atomic_mass)) for sp in species)
            d[file]['radii'] = radii
            d[file]['masses'] = masses
            for sp in species:
         
                if sp.ionic_radius == None:
                # These are charge decorated files which have an assigned oxidation state but no radius in pymatgen corresponding to that state
                # These files will have to be analyzed later, possibly using the crystal radius from the shannon table
                    no_rad.append(file)
                    break

        except:
            # These are files that cannot be assigned oxidation states - bond valence fails either due to disorder or is unable to find a charge neutral state       
            fails.append(file)
            d[file]['struct'] = s
            d[file]['radii'] = None
            d[file]['masses'] = None

    # This is a list of usable files for zeo++
    ready = list(set(files).difference(set(no_rad)).difference(set(fails)))

    for subd in ['zeo-ready', 'zeo-no-rad', 'zeo-fails']:

        # This creates subdirectories if not already present
        if not os.path.isdir('{0}/{1}'.format(dir, subd)):
            os.makedirs('{0}/{1}'.format(dir, subd))

    if remove_duplicates:
        d_ready = remove_duplicates(ready, separator)
        ready = [d_ready[key] for key in d_ready]
        
        d_no_rad = remove_duplicates(no_rad, separator)
        no_rad = [d_no_rad[key] for key in d_no_rad]

        d_fails = remove_duplicates(fails, separator)
        fails = [d_fails[key] for key in d_fails]


    
    # Writing files into respective sub-directories
    for file in ready:
        
        mg.write_structure(d[file]['struct'], '{0}/zeo-ready/{1}'.format(dir, file))
        # write the radius and the mass files also
        filename = file.rsplit('.cif',1)[0]
        write_rad_file(d[file]['radii'],'{0}/zeo-ready'.format(dir), filename)
        write_mass_file(d[file]['masses'],'{0}/zeo-ready'.format(dir), filename) 

    for file in no_rad:
        # These files do not have a radius for atleast one of the species
        mg.write_structure(d[file]['struct'], '{0}/zeo-no-rad/{1}'.format(dir,file))

    for file in fails:
        
        # Note: these files are not charge decorated and are simply the original cif files rewritten
        mg.write_structure(d[file]['struct'], '{0}/zeo-fails/{1}'.format(dir,file))

    return len(ready), len(no_rad), len(fails)

def perform_channel_analysis(dir, prepared = False):

    '''
    Takes a directory and prepares it for analysis by zeo++ (if prepared = False). Then performs channel analysis on all the files in the directory.
    '''
    if not dir.endswith('/'):
        dir+='/'
    if not prepared:
        prepare_for_zeo(dir)
        dir+='ready/'
    files = [file for file in os.listdir(dir) if file.endswith('cif')] 
    structs, sizes = [], []
    for file in files:
        structs.append(file.strip('.cif'))
        sizes.append(find_channel_size(dir + file))
        
    return structs, sizes

   
def cif2cssr(cif, outfile = None, remove = ['Li+']):
    '''
    Converts cif files to Zeo++ CSSR files, deletes species specified in remove. 
    Must have pymatgen installed. Structure must be ordered and oxidation state decorated
    '''
    filename  = cif.rsplit('.',1)[0]
    s = mg.read_structure(cif)

    if remove != None:
        s.remove_species(remove)

    if outfile == None:
        outfile = filename + '.cssr'
        
    try:
        cssr = ZeoCssr(s)
        cssr.write_file(outfile)
    except:
        cssr = None
    return cssr


def find_channel_size(file, accuracy = 'normal'):

    '''
    Use zeo++ to find the largest free sphere. 
    This must be run on a cif file in a directory prepapred for zeo++, i.e., with the radius file and the mass file.
    '''

    filename, ext = file.rsplit('.',1)
    
    if ext == 'cif':
        
        cssr = cif2cssr(file)
        
        if cssr == None:
            
            return None, None
        
        file = '{0}.cssr'.format(filename)
    rad_file = '{0}.rad'.format(filename)
    mass_file = '{0}.mass'.format(filename)
    # Make sure to leave a space on every addition to a command
    cmd = '$ZEO '   

    if accuracy == 'high':
        cmd+='-ha '

    cmd+= '-mass {0} -r {1} -res {2}'.format(mass_file, rad_file, file)
    

    p = Popen(cmd, shell =True, stdout =PIPE, stderr= PIPE)
    out, err = p.communicate()
    if err !='':
        raise Exception('\n\n{0}'.format(err))

    f = open('{0}.res'.format(filename), 'r')   
    free_sp_rad = float(f.readline().split()[2])
    f.close()
    return free_sp_rad


def find_channel_dimensionality(file, probe_radius = 0.5, accuracy = 'normal', use_rad_mass = False):

    '''
    Use zeo++ to find conducting channels in given structure- cif or cssr. Must specify zeo executable as $ZEO in .bashrc for this to work.
    '''

    filename, ext = file.rsplit('.',1)

    if ext == 'cif':

        cssr = cif2cssr(file)

        if cssr == None:

            return None, None

        file = '{0}.cssr'.format(filename)

    cmd = '$ZEO -chan {0}'.format(probe_radius)
    
    if accuracy =='high':

        cmd+=' -ha'
        
    if use_rad_mass:

        cmd+= ' -mass {0}.mass -r {0}.rad'.format(filename)
        
    cmd+= ' {0}'.format(file)

    p = Popen(cmd, shell =True, stdout =PIPE, stderr= PIPE)
    out, err = p.communicate()
    if err !='':
        raise Exception('\n\n{0}'.format(err))

    f = open('{0}.chan'.format(filename), 'r')
    lines = f.readlines()
    f.close()
    channels =  lines[0].split()[1]
    dimensionality = lines[0].split()[-1]
    if int(channels) == 0:
        dimensionality = 0

    return int(channels), int(dimensionality)

    

def write_rad_file(d, path, filename):
    '''
    d = dict of element and ionic radius (element should ideally be charge decorated, but this is not necessary)
    path = path to the directory to store file 
    filename = <filename.rad>
    '''
    
    file = ''
    for key in d:
        file+='{0} {1}\n'.format(key, d[key])

    f = open('{0}/{1}.rad'.format(path, filename), 'w')
    f.write(file)
    f.close
    return


def write_mass_file(d, path, filename):
    '''
    d = dict of element and mass (element should ideally be charge decorated, but this is not necessary)
    path = path to directory to store the file
    filename = file will be stored as <filename.mass>
    '''
    file = ''
    for key in d:
        file+='{0} {1}\n'.format(key, d[key])

    f = open('{0}/{1}.mass'.format(path, filename), 'w')
    f.write(file)
    f.close
    return
    


def read_msd(filepath, skiprows = 1):
    '''
    This function reads a typical msd output file and returns t, msd
    Args : filepath - name of the msd file with two columns t and msd
           skiprows - number of rows to skip on top of the msd file
    Returns: t, msd
    '''

    t, msd = np.loadtxt(filepath, skiprows = skiprows, unpack = True)
    return t, msd



def plot_msd(t, msd,  f = 0.25, skiprows = 0, save = False, show = False, t_units = 'ps', msd_units = '$\AA^{2}$',  slope = True, T = None, legend = False, label = None):

    # TODO make this simpler!!!
    '''
    Plots a mean square displacement trajectory
    Args: save = filepath to save
          show = Turns off/on plt.show()
          Note: Both savefig and show can be done in the main function also
          slope = Turns slope on/off
          T = Temperature for legend
          f = fraction of data to discard for equilibration
          Units are only for the labels
    '''
    
    if label!= None:
        plt.plot(t, msd, label = label)
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

    if legend ==True:
        plt.legend(loc = 'best')

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
        Dint =  np.array(interval)/6. * 1e-4
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


def plot_logD_v_Tinv(Ds, Ts, save = False):


    ln_D = np.log10(np.array(Ds))
    T_inv = 1000./np.array(Ts)

    plt.plot(T_inv, ln_D, 'ro')

    T_inv_stack = np.column_stack([T_inv**1, T_inv**0])

    p, pint, se = regress(T_inv_stack, ln_D, 0.05)

    ln_D_fit = np.dot(T_inv_stack, p)

    plt.plot(T_inv, ln_D_fit)
    plt.xlabel('1000/T (1/K)')
    plt.ylabel('ln(D)')

    if save != False:

        plt.savefig(save)

    slope = p[0]

    R = 5.189e19
    Na = 6.023e23
    E_act = -slope*R/Na
    E_act_int = - np.array(pint[0]) * R/Na

    return E_act, E_act_int


def plot_lnD_v_Tinv(Ds, Ts, save = False):


    ln_D = np.log(np.array(Ds))
    T_inv = 1000./np.array(Ts)

    plt.plot(T_inv, ln_D, 'ro')

    T_inv_stack = np.column_stack([T_inv**1, T_inv**0])

    p, pint, se = regress(T_inv_stack, ln_D, 0.05)

    ln_D_fit = np.dot(T_inv_stack, p)

    plt.plot(T_inv, ln_D_fit)
    plt.xlabel('1000/T (1/K)')
    plt.ylabel('ln(D)')

    if save != False:

        plt.savefig(save)

    slope = p[0]

    R = 5.189e19
    Na = 6.023e23
    E_act = -slope*R/Na

    E_act_int = - np.array(pint[0]) * R/Na


    return E_act, E_act_int

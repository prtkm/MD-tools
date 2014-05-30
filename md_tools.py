
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


def prepare_for_zeo(dir):

    '''
    Takes a directory of cif files and filters zeo++ usable files. Requires pymatgen.
    Creates new directories:
    ./ready - oxidation state decorated files for zeo++ to use.  
    ./no-rad - structures with species having no corresponding ionic radii in pymatgen database
    ./fails - structures which cannot be assigned oxidation states
    '''

    files = [file for file in os.listdir(dir) if file.endswith('cif')]
    
 
    d = {} # Dictionary of the form d[<filename>]['struct'], d[<filename>]['properties']
           
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

    for subd in ['ready', 'no-rad', 'fails']:

        # This creates subdirectories if not already present
        if not os.path.isdir('{0}/{1}'.format(dir, subd)):
            os.makedirs('{0}/{1}'.format(dir, subd))

    # Writing files into respective sub - directories
    for file in ready:
        mg.write_structure(d[file]['struct'], '{0}/ready/{1}'.format(dir, file))
        # write the radius and the mass files also
        filename = file.rsplit('.cif',1)[0]
        write_rad_file(d[file]['radii'],'{0}/ready'.format(dir), filename)
        write_mass_file(d[file]['masses'],'{0}/ready'.format(dir), filename) 

    for file in no_rad:
        mg.write_structure(d[file]['struct'], '{0}/no-rad/{1}'.format(dir,file))

    for file in fails:
        # Note that these files are not charge decorated and are simply the original cif files rewritten
        mg.write_structure(d[file]['struct'], '{0}/fails/{1}'.format(dir,file))

    return

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

def cif2cssr(cif, remove = ['Li+']):
    '''
    Converts cif files to Zeo++ CSSR files, deletes species specified in remove. 
    Must have pymatgen installed. Structure must be ordered and oxidation state decorated
    '''
    filename  = cif.rsplit('.',1)[0]
    s = mg.read_structure(cif)

    if remove != None:
        s.remove_species(remove)
    
    try:
        cssr = ZeoCssr(s)
        cssr.write_file('{0}.cssr'.format(filename))
    except:
        cssr = None
    return cssr


def find_channel_size(file, accuracy = 'normal'):

    '''
    Use zeo++ to find the largest free sphere. 
    This must be run on a cif file in a directory prepapred for zeo++, i.e., with the radius file and the mass file.
    '''

    filename, ext = file.rsplit('.cif',1)
    
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


def find_channels(file, probe_radius = 0.5, accuracy = 'normal', rad_file = None):

# TODO this needs to be rewritten for mass input
    '''
    Use zeo++ to find conducting channels in given structure- cif or cssr. Must specify zeo executable as $ZEO in .bashrc for this to work.
    '''
    filename, ext = file.rsplit('.cif',1)
    
    if ext == 'cif':
        
        cssr = cif2cssr(file)
        
        if cssr == None:
            
            return None, None
        
        file = '{0}.cssr'.format(filename)
      
          
    if accuracy !='high':
        if rad_file == None:           
            cmd = '$ZEO -chan {0} {1}'.format(probe_radius, file)
        else:
            cmd = '$ZEO -nomass -r {2} -chan {0} {1}'.format(probe_radius, file, rad_file)
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

def get_conductivity(atoms, T,  D = None, slope = None, interval = None,  species = 'all'):

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


def plot_lnD_v_Tinv(Ds, Ts, save = False):


    lnD = np.log(np.array(D))
    Tinv = 1000./np.array(Ts)

    plt.plot(lnD, Tinv, 'ro')

    T_inv_stack = np.column_stack([T_inv**1, T_inv**0])

    T_inv_stack = np.column_stack([T_inv_cut**1, T_inv_cut**0])
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

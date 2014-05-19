from ase import *
from numpy import *
import pickle, os, sys, time
from matplotlib import pylab
from pylab import *
from scipy.spatial import Delaunay
from operator import itemgetter
#import psyco
#psyco.full()

'''def prepare_garnet_prototype():
    garnet = Atoms(pbc=True)
    os.chdir("/home/aks1pal/Abinitio/Sneha/Garnet/structure/Analysis")
    for line in open("garnet_prototype.txt"):
        data = line.split()
        garnet += Atom(str(data[0]), ([float(data[1]), float(data[2]), float(data[3])]))
    a = 6.403270
    cell=[[-a,   a,    a],
          [a,   -a,    a],
          [a,    a,   -a]]
    garnet.set_cell(cell, scale_atoms = True)

    Natoms = range(len(garnet))
    Li_atoms = [i for i in Natoms if garnet.get_chemical_symbols()[i]=='Li']
    O_atoms = [i for i in Natoms if garnet.get_chemical_symbols()[i]=='O']
    O_nbrs=[]; O_coord=[]; Li_nbrs=[]
    garnet.site_ids = Li_atoms
    radius = 3.0
    for iLi in Li_atoms:       
        O_near = [iO for iO in O_atoms if garnet.get_min_distance(iLi, iO, mic=True) < 3.0] 
        O_nbrs += [O_near]
        O_coord += [len(O_near)]
#        for jLi in Li_atoms:
        Li_near = set([jLi for jLi in Li_atoms if 0.1 < garnet.get_min_distance(iLi, jLi, mic=True) < 2.5])
        Li_nbrs += [Li_near]
    garnet.Li_nbrs = Li_nbrs
    garnet.O_nbrs = O_nbrs
    garnet.O_coord = O_coord
    write_pickle('garnet_prototype.pkl', garnet)
'''    

def write_pickle(filename, object):
    # pickle the calculation instance for record
    f = open(filename, 'wb')
    pickle.dump(object, f)
    f.close()

    
def read_pickle(filename):
    f = open(filename, 'rb')
    return pickle.load(f)
    f.close()


#prototype = read_pickle("C:\Documents and Settings\kyb5pal\My Documents\Battery\Ion conductors\Structures\MD\garnet_prototype.pkl")
# prototype = read_pickle("garnet_prototype.pkl")
    
def mean_square_displacement(traj):
    # Compute how far each atom moves in an MD run
    Images = range(len(traj))
    Atoms = range(len(traj[0]))
 #  Li_atoms = [i for i in Atoms if traj[0].get_chemical_symbols()[i]=='Li']
 #  dist = [[linalg.norm(traj[t].get_positions()[i] - traj[0].get_positions()[i]) for i in Li_atoms] for t in Images]
    dist = zeros((len(traj), len(traj[0])))
    for t in Images:
        for i in Atoms:
            dist[t,i] = linalg.norm(traj[t].get_positions()[i] - traj[0].get_positions()[i])
    return dist


#This is outdated
def coordination(traj, atom_indices):
    # Compute how many oxygen nearest neighbors there are for each Li
    Nimages = range(len(traj))
    Natoms = range(len(traj[0]))
    Li_atoms = [i for i in Natoms if traj[0].get_chemical_symbols()[i] in ['Li','Al']]
    O_atoms = [i for i in Natoms if traj[0].get_chemical_symbols()[i]=='O']
    Li_coord = zeros((len(Nimages), len(Li_atoms)))
    O_coord = zeros((len(Nimages), len(Li_atoms)))
    for t in Nimages:
        for iLi in atom_indices:
            O_coord[t,iLi] = 0
            Li_coord[t,iLi] = 0
            for iO in O_atoms:
                O_dist = traj[t].get_distance(iLi, iO, mic=True)                
                if O_dist < 3.3 : O_coord[t,iLi] += 1
            if O_coord[t,iLi] == 5: O_coord[t,iLi] = 4
            if O_coord[t,iLi] in [7,8]: O_coord[t,iLi] = 6
        for iLi in atom_indices:
            for jLi in Li_atoms:
                if (O_coord[t,iLi], O_coord[t,jLi]) in [(4,6),(6,4)]:
                    Li_dist = traj[t].get_distance(iLi, jLi, mic=True)
                    if 0 < Li_dist < 3.0: Li_coord[t,iLi] += 1
    return (O_coord, Li_coord)

#This is outdated
def coord_from_matching(traj):
    # Compute how many oxygen nearest neighbors there are for each Li
    start = time.clock()   # time this routine
    radius = 3.2
    Natoms = range(len(traj[0]))
    Li_atoms = [i for i in Natoms if traj[0].get_chemical_symbols()[i] in ['Li','Al']]
    O_offset = 44 - len(Li_atoms) # O atom indices are corrected because of fewer Li ions than prototype
    O_atoms = [i for i in Natoms if traj[0].get_chemical_symbols()[i]=='O']
    for t, image in enumerate(traj):
        print "... Processing image  %s" % t
        O_nbrs=[]; O_coord_raw=[]; O_coord_lim=[]; site_ids=[]; O_coords=[] 
        for iLi in Li_atoms:
            # distance matching to prototype sites
            O_near = [] 
            for iO in O_atoms:
                O_dist = image.get_distance(iLi, iO, mic=True)                
                if O_dist < radius : O_near += [iO + O_offset]
            O_nbrs += [O_near]
            #O_coord_raw += [len(O_near)]
            #if len(O_near) in [3,4,5]:   O_coord_lim += [4]
            #if len(O_near) in [6,7,8]:   O_coord_lim += [6]
            # now assign Li sites to prototype positions
            site_id = -1; O_coord = 0
            p_site_list = list(prototype.site_ids) #copy
            for p_site in p_site_list:   # each site has a list of coordinating O atoms
                #Assign tetra sites only if all tetra connections match, but not more than one extra
                p_nbrhood = prototype.O_nbrs[p_site]
                if len(set(O_near) & set(p_nbrhood)) >= 4:     # if 4 or more nearby sites match a prototype neighborhood
                    if len(p_nbrhood) == 4 and len(O_near) <= 5:    # if matches tetrahedral and not too many other nbrs
                        site_id = p_site
                        O_coord = 4
                        p_site_list.remove(p_site)
                        break
                    elif len(p_nbrhood) == 6 and len(O_near) > 5:    # if matches octahedral 
                        site_id = p_site 
                        O_coord = 6
                        p_site_list.remove(p_site)
                        break
            if site_id < 0 : 
                print "error assigning Li # %s" % iLi   # can use else clause here
            site_ids += [site_id]
            O_coords += [O_coord]
        # ----- end for iLi loop
        image.O_nbrs = O_nbrs
        image.O_coord_raw = O_coord_raw
        image.O_coord_lim = O_coord_lim
        image.site_ids = site_ids
        image.O_coords = O_coords
        image.Li_atoms = Li_atoms
    print "Coordination computations took %s sec" % (time.clock() - start)
    write_pickle(traj[0].filename + ".match.pkl", traj)
    return traj


def coord_from_poly(traj):
    # Compute how many oxygen nearest neighbors there are for each Li
    start = time.clock()   # time this routine
    Natoms = range(len(traj[0]))
    Li_atoms = [i for i in Natoms if traj[0].get_chemical_symbols()[i] in ['Li','Al']]
    O_offset = 44 - len(Li_atoms) # O atom indices are corrected because of fewer Li ions than prototype
    #O_atoms = [i for i in Natoms if traj[0].get_chemical_symbols()[i]=='O']
    for t, image in enumerate(traj):
        print "... Processing image  %s" % t
        O_nbrs=[]; O_coords=[]; site_ids=[]; polyhedra=[]; poly_origins=[]; poly_hulls=[]; Li_nbrs=[]; Li_near=0
        #Li_near = set([jLi for jLi in Li_atoms if 0.1 < image.get_min_distance(iLi, jLi, mic=True) < 2.5]
        #Li_coords[t,iLi]=0
        #precompute shifted polyhedra
        for id in prototype.site_ids:
            nbrhood = prototype.O_nbrs[id]
            # shifting origin of each polyhedron to the first listed O vertex, to get smallest polyhedron
            poly_origins += [nbrhood[0] - O_offset]
            polyhedron = [image.get_min_vector(nbrhood[0] - O_offset, pO - O_offset, mic = True) for pO in nbrhood]
            poly_hull = set(map(frozenset, Delaunay(polyhedron).convex_hull))
            #if len(poly_hull) not in [4,8]:
            #    print "Strange polyhedron of site %s at time %s" % (id, t)
            polyhedra += [polyhedron]
            poly_hulls += [poly_hull]
        
        # compute topological polyhedral coordination for each Li      
        for iLi in Li_atoms:
            site_id = -1; O_coord = 0; nbrhood =[];
            p_site_list = list(prototype.site_ids) #copy
            for p_site in p_site_list:   # each site has a list of coordinating O atoms
                polyhedron = polyhedra[p_site]
                poly_hull = poly_hulls[p_site]
                origin = poly_origins[p_site]
                vertices = polyhedron + [image.get_min_vector(origin, iLi, mic=True)]             
                vert_hull = set(map(frozenset, Delaunay(vertices).convex_hull))
                #if iLi in [1,6]:
                #    print iLi, p_site, polyhedron, vertices
                if poly_hull == vert_hull:     # if true, then iLi is inside the polyhedron!
                    #if site_id > 0 : print "double assigning Li site %s as %s then %s " % (iLi, site_id, p_site)
                    site_id = p_site
                    O_coord = len(polyhedron)
                    nbrhood = prototype.O_nbrs[p_site]
                    p_site_list.remove(p_site)
                    break  # don't allow duplicate site assignments
            if site_id < 0 : 
                print "error assigning Li # %s" % iLi   # can use else clause here
                 #print site_ids
            site_ids += [site_id]
            O_coords += [O_coord]
            O_nbrs += [nbrhood]
            #Li_nbrs += [Li_near]              
        # ----- end for iLi loop
        image.O_nbrs = O_nbrs
        image.site_ids = site_ids
        image.O_coords = O_coords
       #image.Li_nbrs = Li_nbrs                  
    print "Coordination computations took %s sec" % (time.clock() - start)
    #write_pickle(traj[0].filename + ".poly.pkl", traj)
    return traj


def find_occupancy(traj):
    tet_occ = 0
    oct_occ = 0
    for t,image in enumerate(traj):
        tet_occ += sum([1.0 for x in image.O_coords if x==4])
        oct_occ += sum([1.0 for x in image.O_coords if x==6])
    tet_occ = tet_occ/t
    oct_occ = oct_occ/t
    return (tet_occ, oct_occ)   


def find_hops(traj):
    # hop format: (timestep, actual Li index, (id_from, id_to), (coord_from, coord_to))
    # traj should be supplies already with site_ids information
    # initiate the hops with the trivial info from first position
    hops=[(0,iLi,(traj[0].site_ids[iLi], traj[0].site_ids[iLi])) for iLi in range(len(traj[0].site_ids))]
    hops=[]
    image_prev = traj[0]
    for t,image in enumerate(traj):
        comparison = zip(image_prev.site_ids, image.site_ids)
        mismatch = [i for i in range(len(comparison)) if comparison[i][0] != comparison[i][1]]
        if len(mismatch) > 0:  #usually there is only one hop per time step
            hops += [(t, iLi, comparison[iLi]) for iLi in mismatch]
        image_prev = image
    return hops


def find_hops2(traj):
    # hop format: (timestep, actual Li index, (id_from, id_to), (coord_from, coord_to))
    hops=[[] for i in range(len(traj[0]))]
    image_prev = traj[0]
    for t,image in enumerate(traj):
        comparison = zip(image.site_ids, image_prev.site_ids)
        mismatch = [i for i in range(len(comparison)) if comparison[i][0] != comparison[i][1]]
        if len(mismatch) > 0:            
            hops += [(t, id, comparison[id], (image_prev.O_coords[id], image.O_coords[id])) for id in mismatch]           
        image_prev = image
    return hops


def select_suspects(traj, atoms_indices, filename):
    for iatom in atoms_indices:
        traj[0][iatom].symbol = 'C'
        traj[-1][iatom].symbol = 'C'
    write(filename + '_init.cif', traj[0])
    write(filename + '_fin.cif', traj[-1])


def plot_dist_hist(traj, atom_indices):
    # Plot histogram of each Li atom movement
    dist = mean_square_displacement(traj)
    time = range(len(traj))
    p = subplot(111)
    xlabel("Li index", fontsize = 18)
    ylabel("Displacement ($\AA$)", fontsize = 18)
    for t in time:
        p.scatter(atom_indices, dist[t,atom_indices])
    show()


def plot_dist_time(traj, atom_indices):
    #Plot Li displacements versus time
    dist = mean_square_displacement(traj)
    time = range(len(traj))
    p = subplot(111)
    xlabel("time_step", fontsize = 18)
    ylabel("Displacement ($\AA$)", fontsize = 18)
    for iatom in atom_indices:
        p.scatter(time, dist[:, iatom])
    show()


def plot_coord_hist(traj, atom_indices, time):
    #Plot coordinations of Li in versus index
    (O_coord, Li_coord) = coordination(traj, atom_indices)
    time = range(len(traj))
    #p1 = subplot(211)
    #xlabel("Li index", fontsize = 16)
    #ylabel("O coordination", fontsize = 16)
    #for t in time:
    #    p1.scatter(atom_indices, O_coords[t,atom_indices])
    
    p2 = subplot(212)
    xlabel("Li index", fontsize = 16)
    ylabel("Li coordination", fontsize = 16)
    for t in time:
        p2.scatter(atom_indices, Li_coord[t,atom_indices])
    show()
    

def plot_coord_time(traj, atom_indices):
    #Plot coordinations of Li in list versus time
    (O_coord, Li_coord) = coordination(traj, atom_indices)
    time = range(len(traj))
    p1 = subplot(211)
    xlabel("time step", fontsize = 16)
    ylabel("O coordination", fontsize = 16)
    for iatom in atom_indices:
        p1.scatter(time, O_coord[:,iatom])
    p2 = subplot(212)
    xlabel("time step", fontsize = 16)
    ylabel("Li coordination", fontsize = 16)
    for iatom in atom_indices:
        p2.scatter(time, Li_coord[:,iatom])
    show()


def plot_hops_time(hops):
    #Plot coordinations of Li in list versus time
    p1 = subplot(111)
    xlabel("time step", fontsize = 16)
    ylabel("site id", fontsize = 16)
    #NLi = len([hop for hop in hops if hop[0]==0])  # The initial configuration length is just the number of Li atoms
    NLi = 26
    print NLi
    maxtime = sorted(hops, key=itemgetter(0))[-1][0] # Sort by time then read the last hop's time
    for iatom in range(NLi):
        myhops = [hop for hop in hops if hop[1]==iatom]
        times = [0]
        if len(myhops) > 0:
            sites = [myhops[0][2][0]]
            for hop in myhops:
                times.append(hop[0]-1)
                times.append(hop[0])
                sites.append(hop[2][0])
                sites.append(hop[2][1])
            times.append(maxtime)
            sites.append(myhops[-1][2][1])
            p1.plot(times, sites,'-')
    #grid(True, which= 'major')
    show()
    
    
#        def get_min_distance(self, a0, a1, mic=True):
#        """ Correct min distance between atoms in non-cubic cells
#        bkoz: this is a simple search that needs to be optimized, can still fail
#        """
#        R = self.arrays['positions']
#        D = R[a1] - R[a0]
#        # Look for the shortest one in surrounded 3x3x3 cells
#        distances = []
#        for i in (-1, 0, 1):
#            for j in (-1, 0, 1):
#                for k in (-1, 0, 1):
#                    distances.append(np.linalg.norm(D + np.dot(self._cell, np.array([i,j,k]))))
#        return min(distances)
#    
#    
#        def get_min_vector(self, a0, a1, mic=True):
#        """ Correct min distance between atoms in non-cubic cells
#        bkoz: this is a simple search that needs to be optimized, can still fail
#        """
#        R = self.arrays['positions']
#        D = R[a1] - R[a0]
#        # Look for the shortest one in surrounded 3x3x3 cells
#        distances = []; vectors = []
#        for i in (-1, 0, 1):
#            for j in (-1, 0, 1):
#                for k in (-1, 0, 1):
#                    vector = D + np.dot(self._cell, np.array([i,j,k]))
#                    vectors.append(vector)
#                    distances.append(np.linalg.norm(vector))
#        return vectors[np.argmin(distances)]
#                       
    
    
def read_qe_output(filename):
    """Import Espresso-4.2 type file.
    Reads unitcells, atom positions, energies, and forces from the output file.
    """
        
    from ase.calculators.singlepoint import SinglePointCalculator
    
    if isinstance(filename, str): #check if filename is a string
        f = open(filename)
    else:        # Assume it's a file-like object
        f = filename
    
    traj = []   # initiate the trajectory
        
    while True:   #TODO: need to rewrite as for loop with continue and break statements
        line = f.readline()
        if not line: break  
        if 'number of atoms/cell' in line:
            natoms = int(line.split()[4])
        if 'lattice parameter' in line:
            a0 = float(line.split()[4])
        if 'crystal axes' in line:
            # inital cell typically looks like   a(1) = ( -6.403270  6.403270  6.403270 )
            # we extract this three times, assuming Angstrom units (ibrav=0)      
            cell = []
            for i in range(3): 
                data = f.readline().split()    # skip ahead and read new line
                cell += [[float(data[3]), float(data[4]), float(data[5])]]
        if 'site n' in line:  # read initial coordinates
            atoms = Atoms(pbc=True) # initiate first structure, may appear twice in the header
            for i in range(natoms):
                data = f.readline().split()  # skip ahead and read new line
                # typically looks like this. Careful with spaces for >100 atoms!!         
                #   12           Li  tau( 12) = (   0.0000000   3.2016350   4.8024525  )
                atoms += Atom(str(data[1]), ([float(data[6]), float(data[7]), float(data[8])]))
            if 'cryst' in line:
                atoms.set_cell(cell, scale_atoms = True)  # use the latest found cell, in this case the first              
            elif 'a_0' in line:  
                atoms.set_cell(cell, scale_atoms = False) # assume a0=1A
        if 'CELL_PARAMETERS' in line:
            cell = []
            for i in range(3): 
                data = f.readline().split()    # skip ahead and read new line
                cell += [[float(data[0]), float(data[1]), float(data[2])]]
        if 'ATOMIC_POSITIONS' in line:
            atoms = Atoms(pbc=True) # we are starting a new structure
            for i in range(natoms):
                data = f.readline().split()  # skip ahead and read new line
                # typically looks like this:    La       0.250000000   0.875000000   0.625000000          
                atoms += Atom(str(data[0]), ([float(data[1]), float(data[2]), float(data[3])]))
            if 'crystal' in line:
                atoms.set_cell(cell, scale_atoms = True)  # use the latest found cell              
            elif 'angstrom' in line:  
                atoms.set_cell(cell, scale_atoms = False) # assume a0=1A
        if '!    total energy' in line: # read energy after each iteration
            energy = float(line.split()[4])*units.Rydberg
            atoms.set_calculator(SinglePointCalculator(energy,None,None,None,atoms))  # enter energy into the last atoms object
            traj += [atoms]   # update trajectory only now, after all data is gathered, so no positions without energy 
    f.close()
    traj[0].filename = filename
    return traj
# end def read_qe_outp

#this is broken in ASE
#Nimages = 33
#traj = read(filename+".traj", slice(Nimages))

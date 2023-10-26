from __future__ import division
from pymatgen.core import Structure, PeriodicSite, Composition, Element
from pymatgen.transformations.standard_transformations import OxidationStateDecorationTransformation
import sys
try:
    from pymatpro.transformations.advanced_transformations import VoronoiInsertionTransformation
    from pymatpro.transformations.advanced_transformations import OrderingTransformation
except ImportError:
    raise ImportError("Pymatpro not installed on machine")
    sys.exit()


#Sorting of the structure is important to get a predictable arrangment of atoms in structure object
s_orig = Structure.from_file('POSCAR_prim.vasp', sort=True)

oxi = OxidationStateDecorationTransformation({"N":-2, "O":-2, "Fe":3, "F":-1, "H":1})

#4 H per every O
n_h_per_N = 2

'''
Strategy is to create a structure with O as a dummy species replacing N.
We then loop through the structure for each N (i.e., O position) identifying the 4 H sites for each N.
Looping this way generates tetrahedral N-H coordination + doesn't take too much time.
'''

s_orig.replace_species({Element("N"):Element("O")})

'''
VoronoiInsertionTransformation obtains a set of Voronoi sites that are desired at bond length given below
From trial and error, the cut-off radius of 1.05A here gives closest value to experimental N-H bond-length of 0.99A
Note that NH4 tetrahedra are not regular here, likely because of the electrostatics of the structure. Will have to rely on DFT.
Experimental bond length for N-H from: http://hydra.vcp.monash.edu.au/modules/mod2/bondlen.html
'''
vit = VoronoiInsertionTransformation('H1+', n_h_per_N, bond_lengths={'N2-': 0.958}, midpoints=False, interior_only=False)


#Get collection of N (i.e., O sites) in structure
n_indices = []
for i in range(len(s_orig)):
    if s_orig[i].species.almost_equals(Composition("O")):
        n_indices.append(i)


list_of_structures = []
#Now loop through the N-sites one at a time to find nearby H-sites. Note that the H-sites have to be added in cumulatively to get 16 total H
for j in range(len(n_indices)):
    if j == 0:
        #Choose original structure in first iteration
        s_orig.replace(n_indices[j],Composition("N"))
        s_oxi = oxi.apply_transformation(s_orig)
    else:
        #Choose previous structure with added H
        prev_s = list_of_structures[j-1]

        #The way sorting of structure works, H is always added in position before N/O. Hence, adding j*4 to the index value here and adding the check that the species is O indeed.
        temp_index = n_indices[j] + j*n_h_per_N
        if prev_s[temp_index].species.almost_equals(Composition("O")):
            prev_s.replace(temp_index,Composition("N"))
        else:
            print("The structure sorting seems to have gone wrong. The index for loop is not locating the right O to replace with N. Aborting.")
            sys.exit()

        s_oxi = oxi.apply_transformation(prev_s)


    #Apply the Voronoi Insertion Transformation. This yields a (large) number of possible H+ sites.
    s = vit.apply_transformation(s_oxi)


    '''
    Several similar sites are likely to be created by Voronoi insertion.
    So delete duplicate ones that are less than 0.4A away.
    IMPORTANT: Merging sites makes the occupancy of remaining H sites incompatible with what we need. This is fixed in the upcoming for loop.
    '''
    s.merge_sites(0.4, "delete")


    #Collect the set of regular sites (ordered) and H sites (disordered) in the structure to correct the occupancy of all H sites
    sites = []
    h_sites = []

    for site in s:
        if site.is_ordered:
            sites.append(site)
        else:
            h_sites.append(site)


    for s in h_sites:
        #Make sure the individual H sites that are disordered have the right H occupancy
        sites.append(PeriodicSite({'H+': n_h_per_N/len(h_sites)}, s.frac_coords, s.lattice))

    #Create a fresh structure object with right H site occupancies for final ordering
    s = Structure.from_sites(sites)

    '''
    Order the H sites. This uses an Ewald summation to rank structures and yields only the lowest (electrostatic) energy one.
    Similar to OrderDisorder bit written to optimize speed better.
    Should return one structure per for loop step.
    '''
    ot = OrderingTransformation()
    s = ot.apply_transformation(s)


    #Sorting the structure again to ensure the right order of atoms and making sure all H are clubbed into one group in POSCAR.
    s.remove_oxidation_states()
    s = oxi.apply_transformation(s)
    s = s.get_sorted_structure()
    s.remove_oxidation_states()

    list_of_structures.append(s)

    #Write the final structure after adding in all the needed H.
    if j == len(n_indices)-1:
        s.to(filename='POSCAR_H.vasp',fmt='poscar')

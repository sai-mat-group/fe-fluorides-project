from __future__ import division
from pymatgen.core import Structure, PeriodicSite, Composition, Element


#Sorting of the structure is important to get a predictable arrangment of atoms in structure object
s_orig = Structure.from_file('POSCAR.vasp')
prim = Structure.get_primitive_structure(s_orig)
prim.to(filename='POSCAR_prim.vasp',fmt='poscar')


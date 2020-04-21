import numpy as np
import subprocess
import os
import time

################################################
#                                              #
#   Please run this script with the command    #
#                                              #
#   create_geo_relaxed_PBE_D3.py | tee P.xyz   #
#                                              #
################################################

def main():

   # you can change here L to get a larger phosphorene sheet
   L = 4

   a = 3.299
   b = 4.591

   coords_P_phosphorene = np.array([ \
     [ 0.0000000000 ,       2.7006508373 ,       4.1922576756] , 
     [ 0.0000000000 ,       1.8901589613 ,       6.3112478849] ,
     [ 1.6496251748 ,       0.4052459512 ,       6.3112478726] ,
     [ 1.6496251748 ,       4.1855639485 ,       4.1922576688] ])

   x_min = np.min(coords_P_phosphorene[:,0])
   x_max = np.max(coords_P_phosphorene[:,0]) + (L-1)*a
   y_min = np.min(coords_P_phosphorene[:,1])
   y_max = np.max(coords_P_phosphorene[:,1]) + (L-1)*b

   eps = 0.01

   bond_length_H_P = 1.42

   n_H_atoms = 6*L
   n_P_atoms = L**2*4

   natoms = n_H_atoms + n_P_atoms

   with open('P.xyz', 'w') as f:
      print(natoms,"\n", file=f)

   for i_x in range(L):
      for i_y in range(L):

         R_vec = np.array( [ i_x*a, i_y*b, 0] )

         for i_atom in range(np.size(coords_P_phosphorene[:,0])):
             atom_coord = coords_P_phosphorene[i_atom,:] + R_vec
             with open('P.xyz', 'a') as f:
                print("P  ", '{:12.6f}'.format(atom_coord[0]), '{:12.6f}'.format(atom_coord[1]), '{:12.6f}'.format(atom_coord[2]) , file=f)

             if np.abs(atom_coord[0] - x_min) <eps:
                with open('P.xyz', 'a') as f:
                   print("H  ", '{:12.6f}'.format(atom_coord[0]-bond_length_H_P), '{:12.6f}'.format(atom_coord[1]), '{:12.6f}'.format(atom_coord[2]), file=f)
             if np.abs(atom_coord[0] - x_max) <eps:
                 with open('P.xyz', 'a') as f:
                   print("H  ", '{:12.6f}'.format(atom_coord[0]+bond_length_H_P), '{:12.6f}'.format(atom_coord[1]), '{:12.6f}'.format(atom_coord[2]), file=f)
             if np.abs(atom_coord[1] - y_min) <eps:
                 with open('P.xyz', 'a') as f:
                   print("H  ", '{:12.6f}'.format(atom_coord[0]), '{:12.6f}'.format(atom_coord[1]-bond_length_H_P), '{:12.6f}'.format(atom_coord[2]), file=f)
             if np.abs(atom_coord[1] - y_max) <eps:
                 with open('P.xyz', 'a') as f:
                   print("H  ", '{:12.6f}'.format(atom_coord[0]), '{:12.6f}'.format(atom_coord[1]+bond_length_H_P), '{:12.6f}'.format(atom_coord[2]) , file=f)

   os.system("sort -t= P.xyz -o P.xyz")
   os.system("sed -ri '1,2!b;1h;1!H;2!d;x;s/^([^\\n]*)(.*\\n)(.*)/\\3\\2\\1/' P.xyz")

   print("Please fix in the input file all phosphor atoms with LIST", str(n_H_atoms+1)+".."+str(natoms))
   print("Cell size could be: ABC  ", '{:4.1f}'.format(x_max-x_min+15), '{:4.1f}'.format(y_max-y_min+15), 12.0)

if __name__ == "__main__":
    main()

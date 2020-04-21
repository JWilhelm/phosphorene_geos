import numpy as np

def main():

   L = 4

   a = 3.299
   b = 4.591

   coords_P_phosphorene = np.array([ \
     [ 0.0000000000 ,       2.7006508373 ,       4.1922576756] , 
     [ 0.0000000000 ,       1.8901589613 ,       6.3112478849] ,
     [ 1.6496251748 ,       0.4052459512 ,       6.3112478726] ,
     [ 1.6496251748 ,       4.1855639485 ,       4.1922576688] ])

   for i_x in range(L):
      for i_y in range(L):

         R_vec = np.array( [ i_x*a, i_y*b, 0] )

         for i_atom in range(np.size(coords_P_phosphorene[:,0])):
             atom_coord = coords_P_phosphorene[i_atom,:] + R_vec
             print("P  ", '{:12.6f}'.format(atom_coord[0]), '{:12.6f}'.format(atom_coord[1]), '{:12.6f}'.format(atom_coord[2]) )

if __name__ == "__main__":
    main()
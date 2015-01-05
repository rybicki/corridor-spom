import scipy
import grid_model as gm
from model import Species

def compare_matrices(A,B,error):
    assert A.shape == B.shape, "Matrix shape mismatch!"
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            assert A[i,j] == B[i,j], "%s: %i, %i. %s != %s" % (error, i,j, A[i,j], B[i,j])

def random_species(n=5, seed=1234567):
    scipy.random.seed(seed)
    species_list = []
    for i in range(n):
        params = [scipy.random.random() for i in range(5)]
        species_list.append( Species(*params) )
    return species_list


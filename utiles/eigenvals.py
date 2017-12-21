from scipy.linalg import eigh,eig
from pylab import *
def mat_eigh(mat):
	""" Print the eigen  value and vectors for a symmetrix matrix """
	data = matrix(mat)
	eigval, eigvec = eigh(mat)
	print "Eigen values ", eigval
	print "Eigen vectors ", eigval


def mat_eig(mat):
	""" Print the eigen  value and vectors for a matrix """
	data = matrix(mat)
	eigval, eigvec = eig(mat)
	print "Eigen values ", eigval
	print "Eigen vectors ", eigval

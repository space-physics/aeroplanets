import numpy as np

def mat_eigh(mat):
	""" Print the eigen  value and vectors for a symmetrix matrix """
	eigval, eigvec = np.linalg.eigh(mat)
	print("Eigen values ", eigval)
	print("Eigen vectors ", eigval)


def mat_eig(mat):
	""" Print the eigen  value and vectors for a matrix """
	eigval, eigvec = np.linalg.eig(mat)
	print("Eigen values ", eigval)
	print("Eigen vectors ", eigval)

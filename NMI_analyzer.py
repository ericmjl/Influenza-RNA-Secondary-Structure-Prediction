import pandas
import matplotlib.pyplot as plt
from sklearn.metrics import normalized_mutual_info_score as nmi
import numpy as np

def Calculate_NMI(numpy_array):
	"""
	Takes in a numpy array (which has been presumably parsed from an alignment), 
	calculates the normalized mutual information, and returns it as a vector.
	"""
	
	pandas_vector = pandas.DataFrame(numpy_array)
	
	row_length = len(numpy_array[0])
	nmi_vector = np.zeros((row_length,row_length)) # Initialize matrix of zeros.

	# Then fill it with the MI scores.
	for i in range(row_length):
	    for j in range(row_length):
	        if i < j:
	            nmi_vector[i,j] = nmi(pandas_vector.ix[:,i],pandas_vector.ix[:,j]) 
	
	nmi_vector_trans = nmi_vector.T
	np.fill_diagonal(nmi_vector,0)
	nmi_vector_square = nmi_vector + nmi_vector_trans
	np.fill_diagonal(nmi_vector_square,1)
	
	
	return nmi_vector_square
	
	
def Plot_Color(nmi_vector, plot_filename):
	""" 
	Creates a heatmap of the mutual information, prints it to a PDF
	"""
	plt.pcolor(nmi_vector)
	cb = plt.colorbar()
	cb.set_label('Normalized Mutual Information')
	plt.savefig(plot_filename) 
	plt.show()
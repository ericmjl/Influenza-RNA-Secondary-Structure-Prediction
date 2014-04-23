import pandas
import matplotlib.pyplot as plt
from sklearn.metrics import mutual_info_score as mi
import numpy as np
def Calculate_MI(numpy_array):
	"""
	Takes in a numpy array (which has been presumably parsed from an alignment), 
	calculates the mutual information, and returns it as a vector.
	"""
	
	pandas_vector = pandas.DataFrame(numpy_array)
	
	row_length = len(numpy_array[0])
	mi_vector = np.zeros((row_length,row_length)) # Initialize matrix of zeros.

	# Then fill it with the MI scores.
	for i in range(row_length):
	    for j in range(row_length):
	        if i < j:
	            mi_vector[i,j] = mi(pandas_vector.ix[:,i],pandas_vector.ix[:,j]) 
	
	mi_vector_trans = mi_vector.T
	mi_vector_square = mi_vector + mi_vector_trans
	
	
	return mi_vector_square
	
	
def Plot_Binary(mi_vector, plot_filename, cmap='binary'):
	""" 
	Creates a binary heatmap of the mutual information, prints it to a PDF
	"""
	plt.pcolormesh(mi_vector, cmap=cmap)
	# cb = plt.colorbar(orientation='horizontal')
	# cb.set_label('Mutual Information')
	fig = plt.gcf() #get current figure
	fig.set_size_inches(6, 6) #set figure size to 6 inches by 6 inches
	plt.savefig(plot_filename) 
	plt.show()

	return fig
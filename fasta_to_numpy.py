from Bio import AlignIO
import numpy as np

def Parse_and_reverse_complement(fasta_alignment_filename):
	"""
	This function will take a multiple sequence alignment as a fasta file, then parse it using AlignIO, 
	reverse complement it, and return a numpy array with the resulting sequences.
	
	INPUT: 
	fasta_alignment, a multiple sequence alignment as a Fasta file
	
	OUTPUT:
	np_alignment, a reverse-complemented numpy array of the input.
	
	"""
	parsed_alignment = AlignIO.read(fasta_alignment_filename, 'fasta')
	np_alignment = np.array([list(rec.seq.reverse_complement()) for rec in parsed_alignment], np.character)
	
	return np_alignment
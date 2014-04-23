from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

def GetStartingGapLength(sequence):
	"""
	This function takes in a sequence, and counts from the beginning how many 
	gaps exist at the beginning of that sequence, up till the first ATG (start
	codon).

	Parameters:
	- sequence: a string along which you wish to count the number of gaps
	"""
	gaps = 0
	for i, position in enumerate(sequence):
		if sequence[i:i+3] == 'ATG':
			break
		elif sequence[i] == '-':
			gaps += 1
		else:
			pass

	return gaps

def GetEndingGapLength(sequence):
	"""
	This function takes in a sequence, and counts from the end how many 
	gaps exist at the end of that sequence, up till the last letter, which
	could be a stop codon, or not.

	Parameters:
	- sequence: a string along which you wish to count the number of gaps
	"""
	gaps = 0
	for i, position in enumerate(sequence[::-1]):
		if sequence[::-1][i] in IUPACAmbiguousDNA.letters:
			break
		elif sequence[::-1][i] == '-':
			gaps += 1
		else:
			pass

	return gaps

def GetEndGapLengths(sequence):
	"""
	This is syntactic sugar for simultaneously getting the start and end 
	gap lengths.
	"""

	start_gap_length = GetStartingGapLength(sequence)
	end_gap_length = GetEndingGapLength(sequence)

	return start_gap_length, end_gap_length

def GetListOfGapLengths(alignment):
	"""
	This function takes in a multiple sequence alignment and identifies
	the starting gap length for each sequence.

	Parameters:
	- alignment: a BioPython Multiple Sequence Alignment object containing
				 multiple sequences aligned together.
	"""
	gaplengths = []
	for record in alignment:
		sequence = str(record.seq)
		starting_gaps = GetStartingGapLength(sequence)
		ending_gaps = GetEndingGapLength(sequence)
		gaplengths.append((starting_gaps, ending_gaps))

	return gaplengths


def TrimSequence(sequence, front_cutoff, end_cutoff):
	"""
	This function takes a sequence and trims the ends off by the specified
	cutoffs.

	Parameters:
	- front_cutoff: the number of positions to trim off at the front
	- end_cutoff: the number of positions to trim off at the end
	- sequence: the sequence to be trimmed
	"""

	return sequence[front_cutoff:-end_cutoff]

def GetHeaderSequences(sequence, header_positions):
	"""
	This is syntactic sugar for getting the first n nucleotides in a string.

	Parameters:
	- sequence: the sequence of interest
	- header_positions: the number of positions at the head end that you would
						like to get back.
	"""
	return sequence[0:header_positions]

def GetTailSequences(sequence, tail_positions):
	"""
	This is syntactic sugar for getting the last n nucleotides in a string.

	Parameters:
	- sequence: the sequence of interest.
	- tail_positions: the number of positions at the tail end that you would
					  like to get back.
	"""
	return sequence[-tail_positions:]

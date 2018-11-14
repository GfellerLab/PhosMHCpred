#!/usr/bin/env python


######
##
## Function to add calculate peptide score for given PWM and peptide
## Part of PhosMHCpred
##
######


## import packages 
import sys
import math as m


## function to build PWM 
def main(peptide, pwm):

	proteome = { 'A': 7.932417133401348, 
		             'C': 1.5342936224179027, 
		             'E': 6.512641124767595, 
		             'D': 3.481505973643981, 
		             'G': 3.946010725601614, 
		             'F': 4.138140979306901, 
		             'I': 5.323658505848164, 
		             'H': 3.524454881804963, 
		             'K': 6.21126927450847, 
		             'M': 1.7254336310448246, 
		             'L': 10.546829332729441, 
		             'N': 3.130346217109068, 
		             'Q': 4.759669127980483, 
		             'P': 5.494542272076786, 
		             'S': 6.76184763892975, 
		             'R': 6.744164007614165, 
		             'T': 5.247151664757182,
		             'W': 0.6761489469505008, 
		             'V': 7.594570626529927, 
		             'Y': 2.2878896335333394, 
					 's' : 0.5963672127524747, 
					 't' : 0.2533177390191028,
					 'y' : 0.16018072769392427}

	aminoAcids = { 'A' : 0,
					   'C' : 1,
					   'D' : 2,
					   'E' : 3,
					   'F' : 4,
					   'G' : 5,
					   'H' : 6,
					   'I' : 7,
					   'K' : 8,
					   'L' : 9,
					   'M': 10,
					   'N': 11,
					   'P': 12,
					   'Q': 13,
					   'R': 14,
					   'S': 15,
					   'T': 16,
					   'V': 17,
					   'W': 18,
					   'Y': 19,
					   's': 20,
					   't': 21,
					   'y': 22 }
 
	## calculate score
	v = []
	for p in range(len(peptide)):
		v.append( m.log( (pwm[aminoAcids[peptide[p]], p] / float(proteome[peptide[p]]) ) ) )
	score = sum(v) / float( len(v) )
	return score


## call function building PWM 
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])


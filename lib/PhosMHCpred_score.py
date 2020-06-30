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

	proteome = {
	'A' : 6.90102431556,
	'C' : 2.18954997358,
	'E' : 7.03244922217,
	'D' : 4.74227988656,
	'G' : 6.49263488703,
	'F' : 3.54136113332,
	'I' : 4.25444701556,
	'H' : 2.58476567398,
	'K' : 5.63717941422,
	'M' : 2.17252362414,
	'L' : 9.80276546854,
	'N' : 3.52032790345,
	'Q' : 4.73258855788,
	'P' : 6.26735102397,
	'S' : 8.34077088186,
	'R' : 5.60578382853,
	'T' : 5.42923041882,
	'W' : 1.22722764473,
	'V' : 5.94090403552,
	'Y' : 2.58483509058,
	's' : 0.590541123319,
	't' : 0.250843002362,
	'y' : 0.15861587432}


	aminoAcids = { 
	'A' : 0,
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


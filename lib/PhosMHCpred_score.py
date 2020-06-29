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
	'A' : 7.85309296207,
	'C' : 1.51895068619,
	'E' : 6.44751471352,
	'D' : 3.44669091391,
	'G' : 3.90655061835,
	'F' : 4.09675956951,
	'I' : 5.27042192079,
	'H' : 3.48921033299,
	'K' : 6.14915658176,
	'M' : 1.70817929473,
	'L' : 10.4413610394,
	'N' : 3.09904275494,
	'Q' : 4.7120724367,
	'P' : 5.43959684936,
	'S' : 7.28463270317,
	'R' : 6.67672236754,
	'T' : 5.44546470974,
	'W' : 0.669387457481,
	'V' : 7.51862492026,
	'Y' : 3.82656716759,
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


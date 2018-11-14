#!/usr/bin/env python


######
##
##
## PhosMHCpred is a predictor for HLA-I - phosphorylated ligand interactions.
##
## CITATION:
## If you use PhosMHCpred in a publication, please cite:
## Solleder et al., BioRxiv (2018).
##
## FOR-PROFIT USERS:
## If you plan to use PhosMHCpred or any data provided with the script in any for-profit
## application, you are required to obtain a separate  license. To do so, please contact
## eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for  Cancer Research Ltd.
##
## CONTACT:
## Marthe Solleder (marthe.solleder@unil.ch) and David Gfeller (david.gfeller@unil.ch).
##
##
######


## import packages
import argparse, os, sys
from operator import itemgetter
import numpy as np

## import function for predictor
from lib.PhosMHCpred_score import main as scoreF


## path to predictor directory
path = 'YOUR PATH TO /PhosMHCpred FOLDER'
pathTD = path + '/trainingData/'


## inpput arguments
parser = argparse.ArgumentParser()
parser.add_argument( '-p', 
					 '--peptides', 
					 help='Give a list of peptides for prediction. \
					 Accpted format is a text file with one peptide per line or comma-separated in \
					 the command line (i.e. ASEILPPtL,EMDPVtQLY,FTDIEtLKQ)',
					 required=True )
parser.add_argument( '-a', 
					 '--alleles', 
					 help='Give a list of alleles for which the peptides should be predicted, \
					 text file with linewise list of alleles or comma-separated in the command line.',
					 required=True )
parser.add_argument( '-o',
					 '--output',
					 help='Output file. Default: results.txt in current directory.',
					 default='./' )
args = parser.parse_args()

## read peptides
listPeptides = []
if '.txt' in args.peptides:
	listPeptides = open(args.peptides, 'r').read().split('\n')
	listPeptides = list(filter(None, listPeptides))
else:
	listPeptides = args.peptides.split(',')

## checking residues
aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','s','t','y']
for p in listPeptides:
	if len(p) != 9:
		print 'Only 9-mer peptides allowed.'
		print p, 'is length', len(p)
		exit(1)
	for r in p:
		if r not in aa:
			print 'Not supported residue present in peptides: ', r
			print 'Allowed alphabet: ACDEFGHIKLMNPQRSTVWYsty.'
			exit(1)


## length of peptides
length = str(len(listPeptides[0]))

## read / define alleles
listAlleles = []
if '.txt' in args.alleles:
	listAlleles = open(args.alleles, 'r').read().split('\n')
	listAlleles = list(filter(None, listAlleles))
else:
	listAlleles = args.alleles.split(',')
#listAlleles = sorted(listAlleles)

## check if requested alleles are in training data
allelesTD = open( pathTD + 'alleles.txt', 'r' ).read().split('\n')
for x in listAlleles:
	if x not in allelesTD:
		print 'Allele', x, 'not provided in training data...'
		print 'Check alleles.txt in trainingData directory for allowed alleles.'
		exit(1)


## path for output [optional]
if args.output == '':
	pathOutput = './results.txt'
else:
	pathOutput = args.output


## start predictions
print 'Predicting ... '
predictionResults = {}
for hla in listAlleles:
	print hla

	## read pwm	
	pwm = np.loadtxt(pathTD + length + 'mers/' + hla + '_PWM.txt')

	## calculate peptide score
	for peptide in listPeptides:
		if peptide not in predictionResults:
			predictionResults[peptide] = {}
		if hla not in predictionResults[peptide]:
			predictionResults[peptide][hla] = 0
		score = scoreF( peptide, 
						pwm )
		predictionResults[peptide][hla] = score


## print results
output = open(pathOutput, 'w')
output.write('peptide')
for a in listAlleles:
	output.write('\t' + a + '_score')
output.write('\n')
for p in listPeptides:
	output.write(p)
	for a in listAlleles:
		output.write( '\t' + str(predictionResults[p][a]) )
	output.write('\n')
output.close()

print '... done!'


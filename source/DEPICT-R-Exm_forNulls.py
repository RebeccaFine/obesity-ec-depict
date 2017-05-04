### same DEPICT procedure for null exome chip data
### input: list of genes from each permutation in order 
### take same # of genes as was found in the true data

# input list has already been filtered for:
# no repeats, only variants found in height data, only genes seen in DEPICT, only variants seen in nulls
# if other filters desired, create new input list where they have already been applied


import pandas as pd
import numpy as np 
import sys, operator, re, gzip, random, scipy, scipy.stats
from operator import sub
from collections import defaultdict
import ConfigParser
import sys


######################################
######### Read in config file ########
######################################

Config = ConfigParser.ConfigParser()
Config.read(sys.argv[1])

output_dir = Config.get('OUTPUT_SETTINGS','output_dir')
output_label = Config.get('OUTPUT_SETTINGS','output_label')

ordered_permutations_file = Config.get('INPUT_FILES','ordered_permutations_file')

desired_number_genes = int(Config.get('SETTINGS','desired_number_genes'))
recon_gene_sets_file = Config.get('INPUT_FILES','recon_gene_sets_file')

print 'Reading in settings...done'


# read in recon gene sets
recon_gene_sets = pd.read_csv(recon_gene_sets_file, sep = '\t', header = 0, index_col = 0, compression = 'gzip')
#recon_gene_sets = pd.read_csv('../../../data/recon_gene_sets/shortened_recon_gene_sets/3_recon_gene_sets.txt',sep = '\t', header = 0, index_col = 0)

print 'Reading in gene sets...done'

# make name of output file
output_file = output_dir + output_label+ '_' + str(desired_number_genes) + '.nullteststatistic' 

print 'Name of output:', output_file


###########################################################
### calculate null test statistics for each permutation ###
###########################################################

with open(output_file,'w') as output: # open output file
    
    output.write( '\t'.join(recon_gene_sets.columns) + '\n' ) # write out gene set names (in same order they will be iterated through)
    with open(ordered_permutations_file) as nulls: # open ordered permutation file from sortNullsAndOrderGenes)
        for line in nulls: # go through each permutation (one perm per row)
            all_genes = line.strip().split('\t') # split row into gene list
         
            perm_number = all_genes[0] # first item in list is name of permutation
            print 'permutation number:',perm_number
            top_genes = all_genes[1:desired_number_genes + 1] # take the correct # of top genes

            results = [] # store null test statistics [ perm_1, testStat1, testStat2...]
            for geneSet in recon_gene_sets.columns: # go in order of gene sets in recon gene sets
                test_statistic =  recon_gene_sets.loc[top_genes, geneSet].sum() # get test statistic (sum of "top genes" for the given gene set)
                results.append( str(test_statistic) ) # add test statistic to output 

            
            output.write( perm_number + '\t' + '\t'.join(results) + '\n') # write out line to output (perm1 4.5 5.4 4.4...)

          





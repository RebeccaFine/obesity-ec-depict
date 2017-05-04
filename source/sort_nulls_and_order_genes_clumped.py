
"""
Takes in null exome chip data and sorts by p-value.  For each variant, check if it's in the variant annotation file.
Variant annotation file has already been processed for presence in DEPICT, presence in height data, functional consequence,
MAF, etc.
If it is, grab the gene it's in. No repeats allowed.
Output: sorted list of genes for each run.
"""


import pandas as pd
import glob, random, sys
from collections import defaultdict
import ConfigParser

################################
#### Read in config file #######
################################

Config = ConfigParser.ConfigParser()
Config.read(sys.argv[1])

output_dir = Config.get('OUTPUT_SETTINGS','output_dir')
output_label = Config.get('OUTPUT_SETTINGS','output_label')
perm_name = Config.get('OUTPUT_SETTINGS','perm_name')

filtered_variant_annotations = Config.get('INPUT_FILES','filtered_variant_annotations')
null_files = Config.get('INPUT_FILES','null_files')

files_to_read = Config.get('INPUT_FILES','files_to_read') # this specifies which null files to read (for parallelization) -- indexed FROM 1, inclusive at both ends

#SNP_col = Config.get('INPUT_FILES','SNP_col')
#gene_col = Config.get('INPUT_FILES','gene_col')
SNP_col = 'ChromPos'
gene_col = 'Ensembl_ID_in_DEPICT'

# print range of null lines to read
perm_range = map ( int, files_to_read.split('-') )
print 'permutation files to read:', perm_range

perm_start = perm_range[0] - 1
perm_end = perm_range[1] 

# determine names of files to use based on null_file name and perm_start/perm_end
null_stem, null_end = null_files.split('*')

#counters_to_keep = range(perm_start, perm_end) # list of which files to retain (based on given files from files_to_read)

####################################################################
#### Make dictionary with key = variant, value = gene it is in #####
####################################################################

# this is to identify the gene of each variant observed in null 
# also, at this stage, only variants which passed earlier filtering
# will make it through (because the filtered_variants_annotations
# file has already been filtered for nonsyn, presence in current
# EC data/nulls. etc.)

filtered_variant_annotations_df = pd.read_csv(filtered_variant_annotations, sep = '\t', header = 0)  
variant_gene = defaultdict(list) # dict of lists because if a variant has >1 possible annotation, will be stored ( rs1: [geneA], rs2: [geneB,geneC]...}
for (key, val) in filtered_variant_annotations_df[[SNP_col, gene_col]].itertuples(index=False):
   variant_gene[key].append(val)

#with open(filtered_variant_annotations) as complete_annotations:
#    complete_annotations.readline()
#    for line in complete_annotations:
#        name, IlmnSNP, ChromPos, gene, height_MAF, func_region = line.strip().split('\t')
#        print name, gene
#        variant_gene[name].append(gene) #if extra filter desired, apply here

# note for later (will allow specifying of column names instead of going by the number of the column)
#variant_gene = defaultdict(list)
#for (key, val) in df[[SNP_col, gene_col]].itertuples(index=False):
#   variant_gene[key].append(val)

#################################################################################################
#### Open null and sort by p-value; for each ordered variant, write gene (without replacement) ##
#################################################################################################

output_file = output_dir + output_label + '.nullorderedgenes'
print 'output will be written to %s' %output_file

perm_df = [] # list to hold output (will be list of lists, one for each permutation; [ [perm1,ENS1,ENS2,ENS3], [perm2,ENS3,ENS4,ENS5]...])


for i in range(perm_start + 1, perm_end + 1): # go through file numbers specified
    
    null_file = null_stem + str(i) + null_end # create name of file from stem and given number

    null_name = perm_name + '_perm_' + str(i) # add name for null (e.g. pfizer_1)

    null = pd.read_csv(null_file, header = 0, sep = '\s+', engine = 'python') # open each null file

    #null['ChromPos'] = null.loc[:,['CHR','BP']].apply(lambda x: ':'.join(x.dropna().astype(str).values), axis=1) # make Chrom:Pos column

    # remove sex tests
    null_nosex = null[null.TEST=='ADD']

    # sort by p-value
    sorted_null = null_nosex.sort(columns = 'P')

    # make list of ordered SNPs in this null (sorted from best to worst p-value)
    sorted_variants = sorted_null.SNP.tolist()

    sorted_genes = [null_name] # first item in list is name of permutation

    for variant in sorted_variants: # go through each ordered variant
        genes = variant_gene[ variant ] # get gene(s) associated with that variant
        if len(genes) == 1 and genes[0] not in sorted_genes: # if only one gene AND we haven't seen that gene yet before in this list, add to row
            sorted_genes.append( genes[0] )
        elif len(genes) > 1: # if more than one gene, generate random index and take the gene associated with it
            rand_index = random.randint(0, len(variant_gene[variant])-1)
            if genes [ rand_index] not in sorted_genes:
                sorted_genes.append(genes [ rand_index])
                #print variant, genes, rand_index, genes[rand_index]
    
    # add ordered gene list to list of lists for output
    perm_df.append( sorted_genes )
    #counter += 1

#################################
########## write output #########
#################################

output_file = output_dir + output_label + '.nullorderedgenes'

# one line per permutation, each row contains ordered genes
with open(output_file,'w') as output:
    for perm in perm_df:
        output.write( '\t'.join(perm) + '\n')




"""
Compute p-values for gene set enrichment based on list of significant genes from exome chip
Input: single config file
"""

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import numpy as np
import sys, operator, re, gzip, random, scipy, scipy.stats, os
from operator import sub
from collections import defaultdict
import ConfigParser
import matplotlib.pyplot as plt

#### Read in config file

# In[351]:


Config = ConfigParser.ConfigParser()
Config.read(sys.argv[1])
#Config.read('../../runs/height_allMAF_nonSynSplice_EC_Only/height_exome_allMAF_nonSynSplice_DEPICT-R-Exm.cfg')
#Config.read('../../runs/height_allMAF_nonSynSplice_EC_Only_RestrictedGenes/normalizedZScore/height_exome_allMAF_nonSynSplice_DEPICT-R-Exm_RestrictedGenes.cfg')
#Config.read('../../runs/height_MAFLessThan5_nonSynSplice_EC_Only_RestrictedGenes/height_MAFLessThan5_nonSynSplice_EC_Only_RestrictedGenes_DEPICT-R-Exm.cfg')
#Config.read('../../runs/height_allMAF_nonSynSplice_EC_Only/TEST-height_exome_allMAF_nonSynSplice_DEPICT-R-Exm.cfg')

output_dir = Config.get('OUTPUT_SETTINGS','output_dir')
if not output_dir.endswith('/'):
    output_dir = output_dir + '/'
output_label = Config.get('OUTPUT_SETTINGS','output_label')
enrichment_only = Config.getboolean('OUTPUT_SETTINGS','enrichment_only')

include_top_genes = Config.getboolean('OUTPUT_SETTINGS','include_top_genes')
num_genes_output = int( Config.get('OUTPUT_SETTINGS','num_genes_output') )
z_score_cutoff = float( Config.get('OUTPUT_SETTINGS', 'z_score_cutoff') )
gene_style = Config.get('OUTPUT_SETTINGS','gene_style')

if gene_style.lower() == 'hugo':
    gene_style = 'hugo'
elif gene_style.lower() == 'ensembl':
    gene_style = 'ensembl'
elif gene_style.lower() == 'both':
    gene_style = 'both'
else:
    raise Exception('gene_style should say "hugo", "ensembl", or "both"')


if Config.has_option('INPUT_FILES','metacluster_file'):
    metacluster_file = Config.get('INPUT_FILES','metacluster_file') 
    assign_metaclusters = True
else:
    assign_metaclusters = False

recon_gene_sets_file = Config.get('INPUT_FILES','recon_gene_sets_file')
sig_genes_list = Config.get('INPUT_FILES','sig_genes_list')
null_distributions_file = Config.get('INPUT_FILES','null_distributions_file')
num_pval_perm = int( Config.get('INPUT_FILES','num_pval_perm'))
num_FDR_perm = int( Config.get('INPUT_FILES','num_FDR_perm'))
MP_file = Config.get('INPUT_FILES','MP_file')
GO_file = Config.get('INPUT_FILES','GO_file')
PPI_file = Config.get('INPUT_FILES','PPI_file')
ensembl_hugo_file = Config.get('INPUT_FILES', 'ensembl_hugo_file')

# In[352]:

print 'output directory:', output_dir
print 'output label:',output_label
print 'list of significant genes:', sig_genes_list
print 'null distributions:', null_distributions_file
print 'number of pvalue permutations:', num_pval_perm
print 'number of FDR permutations:',num_FDR_perm
print '\n'


########################################################
# Check that all files and directories specified exist #
#######################################################

if os.path.isdir(output_dir) == False:
    raise Exception('Specified output directory does not exist')
if os.path.isfile(recon_gene_sets_file) == False:
    raise Exception('Specified reconstituted gene sets file does not exist')
if os.path.isfile(sig_genes_list) == False:
    raise Exception('Specified list of significant genes does not exist')
if os.path.isfile(null_distributions_file) == False:
    raise Exception('Specified null distributions file does not exist')
if os.path.isfile(MP_file) == False or os.path.isfile(GO_file) == False or os.path.isfile(PPI_file) == False:
    raise Exception('One of the specified gene set name mapping files does not exist')
if os.path.isfile(ensembl_hugo_file) == False:
    raise Exception('The specified Ensembl-HUGO conversion file does not exist')
if os.path.isfile(output_dir + output_label + '_genesetenrichment.txt'):
    raise Exception('A gene set enrichment file with the specified enrichment name already exists...aborting')

if assign_metaclusters:
    if os.path.isfile(metacluster_file) == False:
            raise Exception('The specified metacluster file does not exist')

    if os.path.isfile( output_dir + output_label + '_genesetenrichment_clustersonly.txt'):
        raise Exception('A gene set enrichment clusters file with the specified enrichment name already exists...aborting')
###################################################################
##### Read in reconstituted gene sets, trait exome chip data #####
###################################################################

# read in list of significant genes to analyze
trait_genes_ensembl = []
with open(sig_genes_list) as genefile:
    for line in genefile:
        trait_genes_ensembl.append( line.strip() )


print 'Number of significant genes:',len(trait_genes_ensembl)
print '\n'

# if the number of significant genes read in doesn't show up anywhere in the null file name, throw an exception
if str( len(trait_genes_ensembl) ) not in null_distributions_file:
    raise Exception('Number of significant genes does not match input number used to calculate test statistic (i.e. name of null distribution file does not contain the number of observed significant genes)...aborting')


# read in reconstituted gene sets
recon_gene_sets = pd.read_csv(recon_gene_sets_file, sep = '\t', header = 0, index_col = 0, compression = 'gzip')

# get labels for gene sets
ID_dict = {} # key = gene set name, value = gene set description (if KEGG or REACTOME, same as gene set name)
with open(MP_file) as MP:
    for line in MP:
        ID, description = line.strip().split('\t')[0],line.strip().split('\t')[1]
        ID_dict[ID] = description
with open(GO_file) as GO:
    for line in GO:
        split_line = line.strip().split('\t')
        ID = split_line[0]
        for item in split_line[1:]:
            #print item
            if 'GO:' not in item and item != '':
                description = item
                break
        ID_dict[ID] = description
with open(PPI_file) as PPI:
    for line in PPI:
        ensembl_ID, description = line.strip().split('\t')
        ID_dict[ensembl_ID] = description
        
for gene_set in recon_gene_sets.columns:
    if gene_set not in ID_dict:
        ID_dict[gene_set] = gene_set


##########################################################
############### Calculate test statistic #################
##########################################################

# calculate test statistic for each gene set ( = sum of z-scores in pathway)
geneSetTestStat = {} # key = gene set, value = sum of z-scores (test statistic)
for geneSet in recon_gene_sets.columns:
    test_statistic =  recon_gene_sets.loc[trait_genes_ensembl, geneSet].sum() 
    geneSetTestStat[ geneSet ] = test_statistic

# sort gene sets by increasing test statistic
sorted_geneSetTestStat = sorted(geneSetTestStat.items(), key=operator.itemgetter(1), reverse = True) #sorted list of tuples


###################################################
############### Calculate p-values ################
###################################################


## function for getting p-value by comparing to permutations (used for both actual p-values and FDR calculation):
## get p-values for each gene set: for each z-score sum, subtract off z-score mean from nulls and divide by SD 
## then, convert from adjusted z-score to p-value (based on normal distribution)

def get_pvalue(data_dict, true_data = True): #inputs: dictionary of gene_set:test_stat pairs, True/False for true data versus FDR
        
    if true_data: # if this is for observed data
        geneSetAdjZ = {} # key = adjusted z-score, value = adjusted z-score
        geneSetPVal = {} #key = gene set, value = p-value
        geneSetDirectionality = {} # key = gene set, value = "enriched" or "deenriched"
    else: # if this is for FDR, just need a list of all the p-values
        all_FDR_pval = []
        

    for gene_set, test_statistic in data_dict.iteritems(): #iterate through each gene set
        
        # if observed data, record whether gene set is up- or downregulated
        if true_data:
            if test_statistic < 0:
                geneSetDirectionality[gene_set ] = 'down'
            else:
                geneSetDirectionality[ gene_set ] = 'up'
        
        # get null distribution for the given gene set
        GS_null_distribution = null_distributions_pval[ gene_set ] 
        
        null_mean = np.mean (GS_null_distribution ) # mean of null distribution
        null_sd = np.std ( GS_null_distribution) # sd of null distribution
    
        adjusted_z = (test_statistic - null_mean) / float(null_sd) # from observed data, subtract null mean and null SD to normalize
    
        if true_data:
            geneSetAdjZ[gene_set] = adjusted_z # if for true data, record adjusted z-score
        
        # if looking only at enrichment, only right side of distribution is significant
        if enrichment_only:
            pval = scipy.stats.norm.cdf(-adjusted_z) # to get p-value, take adjusted z-score and get p from normal distribution
            
        # if looking at both enrichment and de-enrichment, take appropriate side of distribution
        else:
            if test_statistic > 0:
                pval = scipy.stats.norm.cdf(-adjusted_z)
            else:
                pval = scipy.stats.norm.cdf(adjusted_z)
           
        # record p-values 
        if true_data: # if actual data, add to dictionary of gene sets and p-values
            geneSetPVal[ gene_set ] = pval
        else: # if FDR, add to list of all p-values
            all_FDR_pval.append(pval)
        
    if true_data:
        return (geneSetAdjZ, geneSetPVal, geneSetDirectionality)
    else:
        return all_FDR_pval


# get nulls for permutation
print 'location of nulls:', null_distributions_file

print 'getting nulls...'
if '.gz' not in null_distributions_file:
    null_distributions = pd.read_csv(null_distributions_file, sep = '\t', header = 0)
else:
    null_distributions = pd.read_csv(null_distributions_file, sep = '\t', header = 0, compression = 'gzip')

null_distributions_pval = null_distributions.iloc[0:num_pval_perm,:] # take as many of the nulls as specificed for p-value calculation
print '...done. %d permutations for p-value calculation' %num_pval_perm


print 'calculating p-values...'


# get adjusted z-scores, p-values, and gene set directionality results 
geneSetAdjZ, geneSetPVal, geneSetDirectionality = get_pvalue(geneSetTestStat, True)

# sort gene sets by adjusted z-score and p-value
geneSetAdjZSorted = sorted(geneSetAdjZ.items(), key=operator.itemgetter(1), reverse = True) # sort by adjusted z-scores
geneSetPValSorted = sorted(geneSetPVal.items(), key=operator.itemgetter(1), reverse = False) # sort by p-value

print '...done'

# make histogram of p-values
plt.hist( geneSetPVal.values(), bins = 30 )
plt.title( 'P-Value Distribution' )
plt.xlabel('P-Value')
plt.ylabel('Count')
figure_name = '%s/%s_ObsPValDistribution.png' %(output_dir, output_label)
plt.tight_layout()
plt.savefig( figure_name )

#################################################################
####################### Calculate FDR ###########################
#################################################################

# take specified number of nulls for FDR calculation
null_distributions_FDR = null_distributions.iloc[-num_FDR_perm:,:] # take from bottom of file


# calculate p-values for specified # of nulls and record all in a list
all_FDR_pval = [] # stores all null p-values
print 'Getting null permutations for FDR...'
for permutation_number, row in null_distributions_FDR.iterrows():
    print permutation_number
    current_FDR_dict = null_distributions_FDR.loc[ permutation_number, : ].to_dict()
    all_FDR_pval += get_pvalue(current_FDR_dict, False)
print '..done\n'
print '\n'

# get list of thresholds to test for FDR

# get list of thresholds to test for FDR

unformatted_thresholds = [] # store thresholds here

min_threshold = 10 ** (int (np.floor( np.log10 ( min (geneSetPVal.values())) ) ) - 1 ) #get order of magnitude of lowest p-value as minimum p-value threshold to test

# add thresholds to list
threshold = min_threshold
unformatted_thresholds.append( threshold )
while threshold < .000001: # up until .0001, add orders of magnitude (10e-20, 10e-21, etc.) 
    threshold = threshold * 10
    unformatted_thresholds.append(threshold)
while threshold < .00001: # up until .0001, add orders of magnitude (10e-20, 10e-21, etc.)
    threshold = threshold + .0000001
    unformatted_thresholds.append(threshold)
while threshold < .0001: # up until .0001, add orders of magnitude (10e-20, 10e-21, etc.)
    threshold = threshold + .000001
    unformatted_thresholds.append(threshold)
while threshold < .001: # until .001, add increments of .00001 
    threshold = threshold + .00001
    unformatted_thresholds.append(threshold) 
while threshold < .01: # until .01, add increments of .0001
    threshold = threshold + .0001
    unformatted_thresholds.append(threshold)
while threshold < .05: # until .05, add increments of .001
    threshold = threshold + .001
    unformatted_thresholds.append(threshold)

thresholds = []
for t in unformatted_thresholds:
    new_t = float('%.3G' %t)
    thresholds.append(new_t)

print 'Number of thresholds to test:', len(thresholds)
print 'Min threshold:', min(thresholds)
print 'Max threshold:', max(thresholds)
print 'Calculating FDRs...\n'

# get FDR q-value associated with each p-value

fdrThresDict = {} # p-value threshold as key, FDR as value

for t in thresholds:
    #print t
    obsCount = sum(x <= t for x in geneSetPVal.values()) # = rank?
    nullCount = sum(x <= t for x in all_FDR_pval)/ float(num_FDR_perm)

    if obsCount == 0: # if observed count is 0, put in 1
        #print t
        fdrThresDict[t] = 1
    
    else:
        fdrThresDict[t] = nullCount/float(obsCount) # to dictionary, add null/observed (= q-value)

    if fdrThresDict[t] > 1: # if it comes out >1, round down to 1
        fdrThresDict[t] = 1


## determine cutoff for FDR <.01, <.05, <.10, <.20 and print

# function for finding highest p-value that meets FDR threshold condition (<.01/<.05/etc.) -- everything with p-value lower than this
# is set to that FDR, even if q-value is a little higher than that (e.g. q = .053, but FDR <.05 if something with a higher p-value 
# has q<.05)
def find_thresh(alpha): 
    FDR_Thresh = -1 # if no thresholds meet condition, set equal to -1 (so everything is higher than this value)
    for threshold in sorted(thresholds): # go through thresholds (which are in order)
        FDR = fdrThresDict[ threshold ] # get FDR for that threshold
        if FDR < alpha: # semi-hack -- assign .01 cutoff to highest value with FDR <.01
            FDR_Thresh = threshold
    return FDR_Thresh


# get FDR thresholds for each alpha
FDR_Thresh_01 = find_thresh(.01) 
FDR_Thresh_05 = find_thresh(.05) 
FDR_Thresh_10 = find_thresh(.10) 
FDR_Thresh_20 = find_thresh(.20)

if FDR_Thresh_01 != -1:
    print 'cutoff for FDR <.01:',FDR_Thresh_01
else:
    print 'no cutoff for FDR <.01'

if FDR_Thresh_05 != -1:
    print 'cutoff for FDR <.05:',FDR_Thresh_05
else:
    print 'no cutoff for FDR <.05'

if FDR_Thresh_10 != -1:
    print 'cutoff for FDR <.10:',FDR_Thresh_10
else:
    print 'no cutoff for FDR <.10'

if FDR_Thresh_20 != -1:
    print 'cutoff for FDR <.20:',FDR_Thresh_20
else:
    print 'no cutoff for FDR <.20'



# for each gene set, find closest threshold to its p-value to approximate its FDR
geneSetFDR = {}
for gene_set, pval in geneSetPVal.iteritems():
    thresholds_more_than_pval = [ t for t in thresholds if t > pval] # make list of thresholds > observed p-value
    if thresholds_more_than_pval != []: # if there exist thresholds > observed p-value, take threshold with minimum distance from observed p-value
        closest = min(thresholds_more_than_pval, key=lambda x:abs(x-pval)) 
        geneSetFDR[ gene_set ] = fdrThresDict [ closest ] # add qvalue to dictionary
    else:
        geneSetFDR[ gene_set ] = 1 # if there are no thresholds > observed p-value, round q-value up to 1
    #print closest
    #break

print '...done\n'

# make histogram of Type I error (based on nulls used for FDR)
plt.hist( all_FDR_pval, bins = 30 )
plt.title( 'Type I Error' )
plt.xlabel('P-Value')
plt.ylabel('Count')
figure_name = '%s/%s_TypeIError.png' %(output_dir, output_label)
plt.tight_layout()
plt.savefig( figure_name )

###################################################################
################ Make output file for results #####################
###################################################################


# append results to list in order of increasing test statistic
results = []
for gene_set, test_statistic in sorted_geneSetTestStat:
    results.append( [gene_set, ID_dict[gene_set], geneSetPVal[ gene_set], geneSetFDR[ gene_set] ,test_statistic, geneSetAdjZ[ gene_set] ] )  
result_df = pd.DataFrame(results)
result_df.columns = ["Original gene set ID", "Original gene set description", "Nominal P value", "qvalue", "Test statistic","Adjusted test statistic"]


# make column for FDR category

def fdr_func(pval): # assign FDRs based on p-value cutoffs calculated earlier, not directly from qvalues
    if pval < FDR_Thresh_01:
        return '<0.01' 
    elif pval < FDR_Thresh_05:
        return '<0.05'
    elif pval < FDR_Thresh_10:
        return '<0.10'
    elif pval < FDR_Thresh_20:
        return '<0.20'
    else:
        return '>=0.20'

result_df['False discovery rate'] = result_df['Nominal P value'].apply(fdr_func)
     

# reorder columns (for Cytoscape purposes)
cols = ["Original gene set ID", "Original gene set description", "Nominal P value", "False discovery rate","qvalue", "Test statistic","Adjusted test statistic"]
result_df = result_df[ cols ]


# sort in order of p-value, then test statistic
sorted_results_df = result_df.sort(columns = ['Nominal P value', 'Adjusted test statistic'], ascending = [1,0]).reset_index(drop=True)

print 'number of gene sets with FDR < .01:', sorted_results_df [ sorted_results_df['Nominal P value'] < FDR_Thresh_01 ].shape[0]
print 'number of gene sets with FDR < .05:', sorted_results_df [ sorted_results_df['Nominal P value'] < FDR_Thresh_05 ].shape[0]
print 'number of gene sets with FDR < .10:', sorted_results_df [ sorted_results_df['Nominal P value'] < FDR_Thresh_10 ].shape[0]
print 'number of gene sets with FDR < .20:', sorted_results_df [ sorted_results_df['Nominal P value'] < FDR_Thresh_20 ].shape[0]
print '\n'

#########################################################################
################## NEW: Get metacluster information #####################
#########################################################################

if assign_metaclusters: 
    metacluster_df = pd.read_csv(metacluster_file, sep = '\t', header = 0)

    # make dictionary with key = original gene set, value = meta-gene set it belongs to
    metacluster_dict = metacluster_df.set_index('Original_ID').MetaGeneSet.to_dict() 

    sorted_results_df['Meta-gene set ID'] = sorted_results_df['Original gene set ID'].map( metacluster_dict )
    sorted_results_df['Meta-gene set description'] = sorted_results_df['Meta-gene set ID'].map( ID_dict )

    # sanity check
    all_gene_sets = sorted_results_df['Original gene set ID'].tolist()
    for gene_set in all_gene_sets:
        if gene_set not in metacluster_dict:
            raise Exception( 'Metacluster file does not contain all gene sets')

#########################################################################
# NEW: add top significant genes and ranks for each gene set to output #
########################################################################

print 'getting top genes for each gene set...'

# make dictionary for Ensembl ID/HUGO ID conversion

ensembl_dict = {} # key = ensembl ID, value = HUGO ID
with open(ensembl_hugo_file) as x:
    for line in x:
        ensembl_id, hugo_id = line.strip().split('\t')
        ensembl_dict[ ensembl_id ] = hugo_id

print 'number of Ensembl IDs with HUGO identifiers:',len(ensembl_dict)
print '\n'



# make smaller version of recon gene sets only containing significant genes

recon_gene_sets_trait_genes_only = recon_gene_sets[ recon_gene_sets.index.isin(trait_genes_ensembl)]
print 'shape of gene sets with significant genes only:',recon_gene_sets_trait_genes_only.shape
print '\n'

# make list of top genes for each gene set
gene_set_top_sig_genes = {} # will contain key: gene set, value = list of significant genes with highest z-score for membership

for gene_set in recon_gene_sets_trait_genes_only.columns: # go through each gene set
    
    # sort gene set by membership z-score
    sorted_geneset = recon_gene_sets_trait_genes_only.loc[:,gene_set].sort_values(ascending = False) 
    
    # get IDs of top gene sets
    ensembl_ids = sorted_geneset[0:num_genes_output].index.tolist() # list of top ensembl ids
    hugo_ids = [ ensembl_dict[gene] for gene in sorted_geneset[0:num_genes_output].index.tolist() ] # list of top HUGO IDs

    # get corresponding z-scores for pathway membership
    top_scores_unrounded = sorted_geneset[0:num_genes_output].tolist()
    
    top_scores = [ round(top_score, 3) for top_score in top_scores_unrounded ] # round to three decimal places


    top_scores_starred = [ str(score) + '*' if score > z_score_cutoff else str(score) for score in top_scores ] # add asterisk for z-scores above specified threshold
    
        
    # get ranks for each gene 
    all_genes = recon_gene_sets[ gene_set].copy()
    all_genes.sort_values(ascending = False, inplace = True) # make sorted copy of recon gene sets with single gene set (all genes, not just significant ones)
    all_genes = pd.DataFrame(all_genes)
    all_genes.loc[:,'Rank'] = [ i + 1 for i in all_genes.reset_index().index.tolist() ] # add column for rank of each gene within the gene set 
    
    top_ranks = []
    for ensembl_id in ensembl_ids:
        top_ranks.append( all_genes.loc [ensembl_id,'Rank'] ) # get rank for each top gene
        
    
    # zip together gene name, z-score, ranks
    if gene_style == 'ensembl': # if style = ensembl, only output ensembl ids
        zipped = zip( ensembl_ids, top_scores_starred, top_ranks ) # list of tuples [(gene set, zscore), (gene set, zscore)..]
    
    elif gene_style == 'hugo': # if style = hugo, hugo ids only
        zipped = zip( hugo_ids, top_scores_starred, top_ranks ) # list of tuples [(gene set, zscore), (gene set, zscore)..]
    
    elif gene_style == 'both': # if style = both, output 'hugo/ensembl'
        combined_ids = []
        for index in range(len(ensembl_ids)):
            #print index
            combined = '%s/%s' %(hugo_ids[index], ensembl_ids[index]) # add slash between hugo and ensembl ID
            combined_ids.append(combined)
        
        zipped = zip(combined_ids, top_scores_starred, top_ranks) # list of tuples [(gene set, zscore, rank), (gene set, zscore, rank)..]
    
    # format to print is geneID:z-score:rank
    formatted = [str(gene) + ' (' + str(zscore) + ',' + str(rank) + ')' for gene,zscore, rank in zipped] 
    
    # add formatted output to dictionary
    gene_set_top_sig_genes[gene_set] = formatted 
    

print 'done'        
print '\n'

print 'Outputting results...congratulations!'

# make data frame containing genes with top z-scores
top_genes_df = pd.DataFrame( gene_set_top_sig_genes ).T
top_genes_df.reset_index(inplace = True)
top_genes_df.rename( columns = {'index':'Original gene set ID'})

headers = ['Original gene set ID'] # make new headers for df
for i in range(num_genes_output):
    header = 'Reconstituted gene set Z score gene %i' % (i+1)
    headers.append(header)

top_genes_df.columns = headers

final_results = pd.merge( sorted_results_df, top_genes_df, on = 'Original gene set ID')

########################################################################################################
#### make version of output with just the metaclusters are listed (with best representative gene set) ##
########################################################################################################

if assign_metaclusters: 
    final_clusters = final_results.groupby('Meta-gene set ID', as_index = False).apply(lambda t: t[t['Nominal P value'] ==t['Nominal P value'].min()]).reset_index(drop = True)

    # sort by p-value
    final_clusters.sort_values(by = ['Nominal P value','Adjusted test statistic'],ascending = [1,0], inplace = True)
    final_clusters.reset_index(inplace = True, drop = True)
    final_clusters.drop_duplicates( ['Nominal P value', 'Meta-gene set ID'] , inplace = True)

    final_clusters.rename( columns = {'Original gene set ID': 'Best representative gene set ID','Original gene set description': 'Best representative gene set description', 'Nominal P value': 'Best representative nominal P value', 'False discovery rate':'Best representative false discovery rate','qvalue':'Best representative qvalue','Test statistic':'Best representative test statistic','Adjusted test statistic':'Best representative adjusted test statistic'}, inplace = True)

    final_clusters_col_order = ['Meta-gene set ID','Meta-gene set description','Best representative gene set ID','Best representative gene set description','Best representative nominal P value','Best representative false discovery rate','Best representative qvalue','Best representative test statistic','Best representative adjusted test statistic'] 


    final_clusters = pd.merge( final_clusters[ final_clusters_col_order], top_genes_df, left_on = 'Best representative gene set ID',right_on = 'Original gene set ID')
    final_clusters.sort_values(by = ['Best representative nominal P value','Best representative adjusted test statistic'], ascending = [1,0], inplace = True)

################################
####### output results ########
###############################

# take care of rounding/formatting
final_results['Nominal P value'] = final_results['Nominal P value'].map(lambda x: '%.4G' % x)
final_results['Test statistic']  = final_results['Test statistic'].round(3)
final_results['Adjusted test statistic']  = final_results['Adjusted test statistic'].round(3)
final_results['qvalue']  = final_results['qvalue'].map(lambda x: '%.3G' % x)

output_name = output_dir + output_label + '_genesetenrichment.txt'
print 'name of output:',output_name 

final_results.to_csv(output_name, sep = '\t', index = False)

if assign_metaclusters:
    
    # take care of rounding/formatting
    final_clusters['Best representative nominal P value'] = final_clusters['Best representative nominal P value'].map(lambda x: '%.4G' % x)
    final_clusters['Best representative test statistic']  = final_clusters['Best representative test statistic'].round(3)
    final_clusters['Best representative adjusted test statistic']  = final_clusters['Best representative adjusted test statistic'].round(3)
    final_clusters['Best representative qvalue']  = final_clusters['Best representative qvalue'].map(lambda x: '%.3G' % x)
    
    cluster_output_name = output_dir + output_label + '_genesetenrichment_clustersonly.txt'
    final_clusters.to_csv(cluster_output_name, sep = '\t', index = False)

print "all finished. You're a rock star!"





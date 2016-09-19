"""
Compute p-values for gene set enrichment based on list of significant genes from exome chip
Input: single config file
"""


import pandas as pd
import numpy as np
import sys, operator, re, gzip, random, scipy, scipy.stats
from operator import sub
from collections import defaultdict
import ConfigParser


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

###################################################################
##### Read in reconstituted gene sets, height exome chip data #####
###################################################################

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


# read in list of significant genes to analyze
height_genes_ensembl = []
with open(sig_genes_list) as genefile:
    for line in genefile:
        height_genes_ensembl.append( line.strip() )


print 'Number of significant genes:',len(height_genes_ensembl)
print '\n'

##########################################################
############### Calculate test statistic #################
##########################################################

# calculate test statistic for each gene set ( = sum of z-scores in pathway)
geneSetTestStat = {} # key = gene set, value = sum of z-scores (test statistic)
for geneSet in recon_gene_sets.columns:
    test_statistic =  recon_gene_sets.loc[height_genes_ensembl, geneSet].sum() 
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
null_distributions = pd.read_csv(null_distributions_file, sep = '\t', header = 0)

null_distributions_pval = null_distributions.iloc[0:num_pval_perm,:] # take as many of the nulls as specificed for p-value calculation
print '...done. %d permutations for p-value calculation' %num_pval_perm


print 'calculating p-values...'


# get adjusted z-scores, p-values, and gene set directionality results 
geneSetAdjZ, geneSetPVal, geneSetDirectionality = get_pvalue(geneSetTestStat, True)

# sort gene sets by adjusted z-score and p-value
geneSetAdjZSorted = sorted(geneSetAdjZ.items(), key=operator.itemgetter(1), reverse = True) # sort by adjusted z-scores
geneSetPValSorted = sorted(geneSetPVal.items(), key=operator.itemgetter(1), reverse = False) # sort by p-value

print '...done'

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

thresholds = [] # store thresholds here

min_threshold = 10 ** (int (np.floor( np.log10 ( min (geneSetPVal.values())) ) ) - 1 ) #get order of magnitude of lowest p-value as minimum p-value threshold to test

# add thresholds to list
threshold = min_threshold
thresholds.append( threshold )
while threshold < .0001: # up until .0001, add orders of magnitude (10e-20, 10e-21, etc.)
    threshold = threshold * 10
    thresholds.append(threshold)
while threshold < .001: # until .001, add increments of .00001 
    threshold = threshold + .00001
    thresholds.append(threshold) 
while threshold < .01: # until .01, add increments of .0001
    threshold = threshold + .0001
    thresholds.append(threshold)
while threshold < .05: # until .05, add increments of .001
    threshold = threshold + .001
    thresholds.append(threshold)


print 'Number of thresholds tested for FDR:', len(thresholds)
print 'Min threshold for FDR:', min(thresholds)
print 'Max threshold for FDR:', max(thresholds)
print '\n'
print 'Calculating FDR for each gene set...'

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

print '...done'

###################################################################
################ Make output file for results #####################
###################################################################


# append results to list in order of increasing test statistic
results = []
for gene_set, test_statistic in sorted_geneSetTestStat:
    results.append( [gene_set, ID_dict[gene_set], geneSetPVal[ gene_set], geneSetFDR[ gene_set] ,test_statistic, geneSetAdjZ[ gene_set] ] )  
result_df = pd.DataFrame(results)
result_df.columns = ["Original gene set ID", "Original gene set description", "Nominal P value", "qvalue", "Test statistic","Adjusted test statistic"]


# make column for FDR category (indicated FDR <.01, <.05, etc.)

def fdr_func(FDR):
    if FDR < .01:
        return '<0.01' 
    elif FDR < .05:
        return '<0.05'
    elif FDR < .10:
        return '<0.10'
    elif FDR < .20:
        return '<0.20'
    else:
        return '>=0.20'

result_df['False discovery rate'] = result_df.qvalue.apply(fdr_func)
    

# reorder columns (for Cytoscape purposes)
cols = ["Original gene set ID", "Original gene set description", "Nominal P value", "False discovery rate","qvalue", "Test statistic","Adjusted test statistic"]
result_df = result_df[ cols ]


# sort in order of p-value, then test statistic
sorted_results_df = result_df.sort(columns = ['Nominal P value', 'Adjusted test statistic'], ascending = [1,0]).reset_index(drop=True)


print 'number of gene sets with FDR < .05:', sorted_results_df [ sorted_results_df.qvalue < .05 ].shape[0]
print 'number of gene sets with FDR < .01:', sorted_results_df[ sorted_results_df.qvalue < .01 ].shape[0] 
print '\n'

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

recon_gene_sets_height_genes_only = recon_gene_sets[ recon_gene_sets.index.isin(height_genes_ensembl)]
print 'shape of gene sets with significant genes only:',recon_gene_sets_height_genes_only.shape
print '\n'

# make list of top genes for each gene set
gene_set_top_sig_genes = {} # will contain key: gene set, value = list of significant genes with highest z-score for membership

for gene_set in recon_gene_sets_height_genes_only.columns: # go through each gene set
    
    # sort gene set by membership z-score
    sorted_geneset = recon_gene_sets_height_genes_only.loc[:,gene_set].sort_values(ascending = False) 
    
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
top_genes_df.loc[:,'index'].rename('Original gene set ID')

headers = ['Original gene set ID'] # make new headers for df
for i in range(num_genes_output):
    header = 'Reconstituted gene set Z score gene %i' % (i+1)
    headers.append(header)

top_genes_df.columns = headers

# merge original results with top z-scores 
final_results = pd.merge( sorted_results_df, top_genes_df, on = 'Original gene set ID')

# take care of rounding/formatting
final_results['Nominal P value'] = final_results['Nominal P value'].map(lambda x: '%.4G' % x)
final_results['Test statistic']  = final_results['Test statistic'].round(3)
final_results['Adjusted test statistic']  = final_results['Adjusted test statistic'].round(3)
final_results['qvalue']  = final_results['qvalue'].map(lambda x: '%.3G' % x)


################################
####### output results ########
###############################

output_name = output_dir + output_label + '_genesetenrichment.txt'
print 'name of output:',output_name 

final_results.to_csv(output_name, sep = '\t', index = False)


print "all finished. You're a rock star!"





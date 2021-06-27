#usage: python Pathsim.py $drugnum
#input: inputfiles were produced by extracRelation.pl, including drug.txt target.txt relation.txt term.txt
#output: result.txt

import pandas as pd
import numpy as np
import sys

star=int(sys.argv[1])
print(star)
########################################################################################################################
# Load the data as a Pandas DataFrame

drugs = pd.read_csv('drug.txt', sep='\t', header= None, names = ['drug_id', 'drug'])
targets  = pd.read_csv('target.txt', sep='\t', header= None, names = ['target_id', 'target'])
relations = pd.read_csv('relation.txt', sep='\t', header= None, names = ['target_id', 'relation_id', 'score'])
terms = pd.read_csv('term.txt', sep='\t', header= None, names = ['term_id', 'term'])

########################################################################################################################
# Create the commuting matrix M(i,:) and perform reordering
# Parameters: adjacency matrices for the meta-path and drug ID
def commute_matrix(adj1, adj2, drug_id):
    M_p = adj1.dot(adj2)

    order = M_p.index.tolist()
    order.insert(0, order.pop(order.index(drug_id)))
    M_pi = M_p.reindex(order)

    M = M_pi.dot(M_p.transpose())
    mid = M[drug_id]
    M.drop(labels=[drug_id], axis=1,inplace = True)
    M.insert(0, drug_id, mid)
    return M

# Get the scaling terms to compute the PathSim scores
# Parameters: commuting matrix and drug id
def scaling(M, drug_id):
    # Diagonal vector
    D = list(np.diagonal(M))
    scale = []
    for i in D:
        scale.append(2 / (D[0] + i))
    CF = pd.DataFrame(M.ix[drug_id])
    CF.columns = ['pathsim_score']
    s = pd.Series(scale)
    return CF, s

# Utility function to perform element-wise multiplication between Dataframes and Series
mult = lambda x: np.asarray(x) * np.asarray(s)

########################################################################################################################
#  Create the drug-target adjacency

merged = relations.merge(drugs, left_on='relation_id', right_on='drug_id')
drug_target = merged[['target_id', 'drug_id']]
drug_target_adj = pd.DataFrame(0, index=drugs['drug_id'], columns=targets['target_id'])
for target, drug in zip(drug_target['target_id'], drug_target['drug_id']):
    drug_target_adj.set_value(drug, target, 1)
drug_target_adj = drug_target_adj.dropna(axis='columns')


# Create the target-Venue adjacency
merged = relations.merge(terms, left_on='relation_id', right_on='term_id')
target_venue = merged[['target_id', 'term_id']]
target_venue_adj = pd.DataFrame(0, index=targets['target_id'], columns=terms['term_id'])
for target, venue in zip(target_venue['target_id'], target_venue['term_id']):
    target_venue_adj.set_value(target, venue, 1)
target_venue_adj = target_venue_adj.dropna(axis='index')

# #  Create the target-Term adjacency
# merged = relations.merge(terms, left_on='relation_id', right_on='term_id')
# target_term = merged[['target_id', 'term_id']]
# target_term_adj = pd.DataFrame(0, index=targets['target_id'], columns=terms['term_id'])
# for target, term in zip(target_term['target_id'], target_term['term_id']):
    # target_term_adj.set_value(target, term, 1)
# target_term_adj = target_term_adj.dropna(axis='index')

########################################################################################################################
# Kristin Faison id: 55154 | use meta-path drug-target-term-target-term
# Master - Forest and Wood Sciences | PhD - English and Speech Teacher Education

# Create output files
#anhai_pathsim = open('anhai_pathsim.txt', 'w')

########################################################################################################################
# Get the scores
M_anhai = commute_matrix(drug_target_adj, target_venue_adj, star)
AD, s = scaling(M_anhai, star)
an_hai = AD.apply(mult)
an_hai_result = drugs.merge(an_hai, left_on='drug_id', right_index =True)
an_hai_result = an_hai_result.sort_values('pathsim_score', ascending=False)
#print(an_hai_result.head(12), file = anhai_pathsim)
an_hai_result.to_csv('result.txt', sep='\t', header= None)
#print(an_hai_result, file = anhai_pathsim)




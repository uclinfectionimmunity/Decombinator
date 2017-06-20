########################################################################################################################
# Short read decombinator
#
########################################################################################################################

from __future__ import division
from Bio import SeqIO
import string
import time
from acora import AcoraBuilder
import math
from Bio import pairwise2
########################################################################################################################


########################################################################################################################
def setup(vregions_file, jregions_file):

    v_end_length = 40  # how many nts at the end of the V region to consider
    j_start_length = 40  # how many nts at the start of the J region to consider

    handle = open(vregions_file, 'r')
    v_list = list(SeqIO.parse(handle, 'fasta'))
    handle.close()
    v_genes = [str(string.upper(v.seq)) for v in v_list]
    v_genes_cut = [v[-v_end_length:] for v in v_genes]

    all_v_substrings = []
    for v in v_genes_cut:
        all_v_substrings.append([v[i:i+n] for n in range(4, len(v)+1) for i in range(len(v)-(n-1))])

    t0 = time.time()
    v_keyword_tries = []
    for v_substrings in all_v_substrings:
        v_builder = AcoraBuilder()
        for i in range(len(v_substrings)):
            v_builder.add(v_substrings[i])
        v_keyword_tries.append(v_builder.build())
    print 'V keyword tries built in', round(time.time() - t0, 2), 'seconds'

    handle = open(jregions_file, 'r')
    j_list = list(SeqIO.parse(handle, 'fasta'))
    handle.close()
    j_genes = [str(string.upper(j.seq)) for j in j_list]
    j_genes_cut = [j[:j_start_length] for j in j_genes]

    all_j_substrings = []
    for j in j_genes_cut:
        all_j_substrings.append([j[i:i+n] for n in range(4, len(j)+1) for i in range(len(j)-(n-1))])

    t0 = time.time()
    j_keyword_tries = []
    for j_substrings in all_j_substrings:
        j_builder = AcoraBuilder()
        for i in range(len(j_substrings)):
            j_builder.add(j_substrings[i])
        j_keyword_tries.append(j_builder.build())
    print 'J keyword tries built in', round(time.time() - t0, 2), 'seconds'

    return v_keyword_tries, j_keyword_tries, v_genes, j_genes
########################################################################################################################


########################################################################################################################
def analyse_file(infile, outpath, fileid, v_keyword_tries, j_keyword_tries, v_genes, j_genes, parameters):

    print 'Classifying sequences in', infile, 'using parameter set', parameters

    outfile = outpath + fileid + '_beta.txt'
    outhandle = open(outfile, 'w')

    t0 = time.time()
    record_count = 0
    assigned_count = 0

    inhandle = open(infile, 'r')
    for record in SeqIO.parse(inhandle, 'fastq'):

        assigned_j = None
        assigned_v = None
        v_deletions = None
        j_deletions = None
        insert_nt = None

        record_count += 1
        if record_count % 1000 == 0:
            print record_count

        # print 'assigning a J gene'
        j_substrings_by_gene = find_all_matches(record, j_keyword_tries)
        assigned_j = assign_gene(record, j_substrings_by_gene, j_genes, parameters)
        # print assigned_j

        if not assigned_j:
            continue  # move to next record in the file

        # print 'assigning a V gene'
        v_substrings_by_gene = find_all_matches(record, v_keyword_tries)
        assigned_v = assign_gene(record, v_substrings_by_gene, v_genes, parameters)
        # print assigned_v

        if not assigned_v:
            continue  # move to the next record in the file

        # print 'finding V deletions'
        v_deletions = find_v_deletions(v_substrings_by_gene, assigned_v, v_genes)
        # print v_deletions

        # print 'finding J deletions'
        j_deletions = find_j_deletions(j_substrings_by_gene, assigned_j, j_genes)
        # print j_deletions

        # print 'finding inserted nucleotides'
        insert_nt = find_insert_nts(record, v_substrings_by_gene, assigned_v, j_substrings_by_gene, assigned_j)
        # print insert_nt

        # print 'writing dcr to outfile'
        dcr = ','.join([str(assigned_v), str(assigned_j), str(v_deletions), str(j_deletions), insert_nt])
        print >> outhandle, dcr

        assigned_count += 1

    inhandle.close()
    outhandle.close()

    print '-----------------------------------------------------------------------'
    print '{0:,}'.format(record_count), 'sequence reads analysed in', round(time.time()-t0, 2), 'seconds'
    print '{0:,}'.format(assigned_count), 'successfully assigned a classifier'
    print '-----------------------------------------------------------------------'


########################################################################################################################


########################################################################################################################
def find_all_matches(record, keyword_tries):

    substrings_found_by_gene = [trie.findall(str(record.seq)) for trie in keyword_tries]

    return substrings_found_by_gene
########################################################################################################################


########################################################################################################################
def assign_gene(record, substrings_found_by_gene, genes, parameters):

    if not substrings_found_by_gene:
        assigned_gene = None
        return assigned_gene
    else:
        # attempt to match by longest substring
        assigned_gene = assign_by_longest_substring(substrings_found_by_gene, parameters)
        if assigned_gene:
            return assigned_gene
        else:
            # attempt to match by scoring
            assigned_gene = assign_by_gene_scores(substrings_found_by_gene, parameters)
            if assigned_gene == 'NA':
                return None
            elif assigned_gene:
                return assigned_gene
            else:
                # do pairwise alignment on two top scoring genes
                assigned_gene = assign_by_pairwise(substrings_found_by_gene, record, genes)
                return assigned_gene
########################################################################################################################


########################################################################################################################
def assign_by_longest_substring(substrings_found_by_gene, parameters):
    """
    Checks whether the gene with the longest substring match is adequate to be assigned.
    Returns matched gene if adequate, None otherwise
    """

    # print 'attempting assignment by longest substring'

    match_length_threshold = parameters[0]
    match_length_differential = parameters[1]

    longest_substring_by_gene = []
    for matches in substrings_found_by_gene:
        lengths = [len(x) for x in matches]
        longest_substring_by_gene.append(max(lengths))

    if max(longest_substring_by_gene) < match_length_threshold:
        return None

    elif sorted(longest_substring_by_gene, reverse=True)[0] \
            < sorted(longest_substring_by_gene, reverse=True)[1] + match_length_differential:
        return None

    else:
        assigned_gene = longest_substring_by_gene.index(max(longest_substring_by_gene))
        return assigned_gene
########################################################################################################################


########################################################################################################################
def assign_by_gene_scores(substrings_found_by_gene, parameters):
    """
    Attempts to match a gene to the sequence read by scoring each gene according to match lengths
    Returns gene index if matched successfully
            None if cannot match successfully
            'NA' if this read should be discarded due to all scores being too low
    """

    # print 'attempting assignment by scoring'

    score_threshold = parameters[2]
    score_differential = parameters[3]

    gene_scores = []
    for gene_matches in substrings_found_by_gene:
        gene_scores.append(sum([math.e**len(x[0]) for x in gene_matches]))

    if max(gene_scores) < score_threshold:
        assigned_gene = 'NA'
        return assigned_gene

    elif sorted(gene_scores, reverse=True)[0] < sorted(gene_scores, reverse=True)[1] * score_differential:
        assigned_gene = None
        return assigned_gene

    else:
        assigned_gene = gene_scores.index(max(gene_scores))
        return assigned_gene
########################################################################################################################


########################################################################################################################
def assign_by_pairwise(substrings_found_by_gene, record, genes):

    # print 'attempting assignment by pairwise alignment'

    gene_scores = []
    for gene_matches in substrings_found_by_gene:
        gene_scores.append(sum([math.e**len(x[0]) for x in gene_matches]))

    gene_indices_to_test = [i for i, x in enumerate(gene_scores) if x in sorted(gene_scores, reverse=True)[:2]]
    pw_scores = []
    for i in gene_indices_to_test:
        pw_scores.append(max([a[2] for a in pairwise2.align.globalxs(str(record.seq)[3:], genes[i], -0.5, -0.5)]))

    assigned_gene = gene_indices_to_test[pw_scores.index(max(pw_scores))]
    return assigned_gene
########################################################################################################################


########################################################################################################################
def find_v_deletions(substrings_found_by_gene, assigned_gene, genes):

    relevant_substrings = substrings_found_by_gene[assigned_gene]

    maximum_matched_length = max([len(x[0]) for x in relevant_substrings])

    relevant_match = [x[0] for x in relevant_substrings if len(x[0]) == maximum_matched_length][0]

    relevant_gene = genes[assigned_gene]

    position = string.find(relevant_gene, relevant_match)

    v_deletions = len(relevant_gene) - (position + len(relevant_match))

    return v_deletions
########################################################################################################################


########################################################################################################################
def find_j_deletions(substrings_found_by_gene, assigned_gene, genes):

    relevant_substrings = substrings_found_by_gene[assigned_gene]

    maximum_matched_length = max([len(x[0]) for x in relevant_substrings])

    relevant_match = [x[0] for x in relevant_substrings if len(x[0]) == maximum_matched_length][0]

    relevant_gene = genes[assigned_gene]

    position = string.find(relevant_gene, relevant_match)

    j_deletions = position

    return j_deletions
#######################################################################################################################


########################################################################################################################
def find_insert_nts(record, v_substrings_by_gene, v_assigned, j_substrings_by_gene, j_assigned):

    v_max_match_length = max([len(x[0]) for x in v_substrings_by_gene[v_assigned]])
    v_longest_match = [x[0] for x in v_substrings_by_gene[v_assigned] if len(x[0]) == v_max_match_length][0]

    j_max_match_length = max([len(x[0]) for x in j_substrings_by_gene[j_assigned]])
    j_longest_match = [x[0] for x in j_substrings_by_gene[j_assigned] if len(x[0]) == j_max_match_length][0]

    v_end_position = string.find(str(record.seq), v_longest_match) + len(v_longest_match)
    j_start_position = string.find(str(record.seq), j_longest_match)

    insert_nts = str(record.seq)[v_end_position: j_start_position]

    return insert_nts
#######################################################################################################################


########################################################################################################################
if __name__ == '__main__':

    vfile = '/home/niclas/Documents/decombinator-v2.0/human_TRBV_region.fasta'
    jfile = '/home/niclas/Documents/decombinator-v2.0/human_TRBJ_region.fasta'

    v_key, j_key, v_regions, j_regions = setup(vfile, jfile)

    infile = 'testseqs.fastq'
    outpath = '/home/niclas/Documents/decombinator-v2.0/'
    fileid = 'MyShortReadAnalysis'
    param_set = [10, 2, 1400, 1.05]
    analyse_file(infile, outpath, fileid, v_key, j_key, v_regions, j_regions, param_set)
########################################################################################################################

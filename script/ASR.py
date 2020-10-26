'''
@sungml92

ASR v.1b ... last updated on May 14, 2020
annotate input tree for condonML based on observed changes along branches in rst file

Upcoming features in next updates:

1) Easier to read codes
2) Full ASR automation by from Bio.Phylo.PAML import codeml
'''

import re
import itertools
from Bio import Phylo
from io import StringIO
import argparse
from time import time

################################################################################
## Ancestral sequence reconstruction (ASR) v.1 beta
## from PAML software output
## Step 4-b in Tool for Bioinformatics Analysis of virus genome (TBA) pipeline
################################################################################
## Module contents
###################
##
## 1. Matching functions
##
##      patterns finding for Node, Substitution, and Branch
##
## 2. rst file reader
##
##      file reading to handle rst output from PAML-codeml
##
## 3. Node-Substitution assigner
##
##      utility to assign node to substitution dictionary
##
## 4. Output file handling
##
##      output file naming
##
## 5. Main function
##
##      main function for command line tools
##
################################################################################

################################################################################
## Matching functions
######################

def get_ancestral_node(text):
    P = "[A-Za-z\s\t]+([0-9]+)[A-Za-z\t\s]+([0-9]+)[a-z\s\t]+"
    T = text
    A = re.search(P,T)
    ans_start = A[1]
    ans_end = A[2]

    return ans_start, ans_end

def get_substitution(line):
    P = "([0-9]+)[\sA-Z]+\(([A-Z\*])\)[\s0-9\.\-\>A-Z\?]+\(([A-Z\*\-])\)"
    T = line
    A = re.search(P,T)
    site, ans_aa, off_aa = A[1], A[2], A[3]
    return site, ans_aa, off_aa

def get_branch_information(line):
    P = "Branch\s+([0-9]+)\:\s+([0-9]+)..([0-9]+)"
    T = line
    A = re.search(P,T)
    branch_number = A[1]
    ans_node = A[2]
    off_node = A[3]

    return branch_number, ans_node, off_node

################################################################################
## rst file reader
###################

def read_asr(file):

    found_branch = False
    curr_branch = 0

    branch_info = {}
    substitution_info = {}

    for line in file:
        line = line.strip()

        ## retrieve ancestral node start and end position

        if line.startswith("Nodes"):
            start, end = get_ancestral_node(line)

        ## retrieve phylogenetic tree with node annotation

        if line.lower().endswith("treeview"):
            tree = next(file)

        ## retrieve branch change
        ## only start to print when Branch is detected, then stop after all branches are parsed

        if line.startswith("Branch"):
            found_branch = True

        elif line.startswith("List"): ## When List is detected stop printing line
            found_branch = False

        if found_branch: ## Now line will be printed by this, write desired script here

            if line.startswith("Branch"):
                branch_no, ans_node, off_node = get_branch_information(line)
                branch_info[str(branch_no)] = [ans_node,off_node]

                curr_branch = branch_no


            else:
                if len(line) > 0:
                    site, ans_aa, off_aa = get_substitution(line)
                    if ans_aa != off_aa:
                        change = ans_aa+site+off_aa
                        if change is not None:
                            try:
                                substitution_info[str(curr_branch)].append(change)
                            except KeyError:
                                substitution_info[str(curr_branch)] = [change]

    return start, end, tree, branch_info, substitution_info

################################################################################
## Node-Substitution assigner
##############################

def gNS(start, branch_info, sub_info):
    '''
    generate Node-Substitution dictionary
    '''

    branch_node_info = {branch: node_tip[1] for branch, node_tip in branch_info.items() if int(node_tip[0]) > int(start)-1 and int(node_tip[1]) > int(start)-1}
    common_keys = list(set(branch_node_info.keys()) & set(sub_info.keys()))
    filtered_bn_info = {branch: node for branch, node in branch_node_info.items() if branch in common_keys}
    filtered_sub_info = {branch: sub for branch, sub in sub_info.items() if branch in common_keys}

    node_sub_info = dict()

    for (branch1, node), (branch2, sub) in zip(filtered_bn_info.items(), filtered_sub_info.items()):
        if branch1 == branch2:
            node_sub_info[node] = sub

    return node_sub_info

################################################################################
## Output file handling
########################

def make_asr_tree(rst_file,reference_tree,mode,out):

    rst_f = open(rst_file)
    ref_t = open(reference_tree,"r").read()

    start, end, tree, branch_info, sub_info = read_asr(rst_f)
    node_sub_dict = gNS(start,branch_info,sub_info)

    ## modify node information in nexus style

    for node in range(int(start), int(end)+1):
        tree = tree.replace(" %s "%node, "[&node-number=%s]"%node)

    ## use replace to fill phylogenetic tree

    for node, sub in node_sub_dict.items():
        sub = '-'.join(sub)
        tree = tree.replace("[&node-number=%s]"%node,"[&node-number=%s,aa-sub=%s]"%(node,sub))

    ## extract branch length from reference_tree (tree used for PAML analysis)

    pos = 0
    b_lengths = []

    while pos < len(ref_t):
        blength = re.match("\:[\-0-9\.E]+",ref_t[pos:pos+50])
        if blength is not None:
            pos += int(blength.span()[1])
            b_lengths.append(blength.group())
        elif blength is None:
            pos += 1

    ## find position to put branch_length

    hanger = re.findall("([0-9]+\_B\/[0-9A-Za-z\-\_\.\|\/]+|\[[\&A-Za-z0-9\-\=\,]+\])",tree)

    if mode == "rooted":
        for hplacement, blength in zip(hanger[:-1],b_lengths):
            tree = tree.replace(hplacement,hplacement+blength)

    elif mode == "unrooted":
        for hplacement, blength in zip(hanger,b_lengths):
            tree = tree.replace(hplacement,hplacement+blength)

    ## clean tip number

    tree = re.sub("[0-9]+\_B/","B/",tree)

    ## read tree as newick and save as nexus

    tree = Phylo.read(StringIO(tree),"newick")

    with open(out.split(".")[0]+".asr.tre", "w") as f:
        Phylo.write(tree,f,"nexus")

################################################################################
## Main function
#################

def main():

    start_time = time()

    parser = argparse.ArgumentParser()
    parser.add_argument("--in_rst","-ir", dest="rst_file", help="input rst file")
    parser.add_argument("--in_tree","-it", dest="reference_tree", help="input reference tree")
    parser.add_argument("--mode","-m", dest="mode", help="reference tree mode, rooted or unrooted")
    parser.add_argument("--out","-o", dest="out", help="output annotated tree")

    args = parser.parse_args()

    make_asr_tree(args.rst_file, args.reference_tree,
        args.mode, args.out)

    end_time = time()

    tot_time = end_time - start_time
    print('\n** Total elapsed runtime:', str(int((tot_time / 3600))) + ':' + str(int((tot_time % 3600) / 60))
         + ':' +str(int((tot_time % 3600) % 60)))

if __name__ == "__main__":
    main()

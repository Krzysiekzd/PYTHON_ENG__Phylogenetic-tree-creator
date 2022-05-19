import datetime
import numpy as np
from typing import List, Tuple
from Bio import SeqIO
import random


class Sample:
    def __init__(self, id: str, country: str, seq: str, date: datetime.date):
        self.sample_id = id
        self.country = country
        self.seq = seq
        self.date = date

'''
Node extends the Sample class. It is a structure belonging to a tree. A node has links to its children in the form of the children attribute,
which is a list of child nodes. If the list is empty, the node is a leaf.
'''
class Node:
    def __init__(self, sample: Sample):
        self.sample = sample
        self.children = []

class Tree:
    def __init__(self, root: Node):
        self.root = root

    # Methods that returns list of all the edges of the tree. Each edge is a pair (id_of_the_parent_sample, id_of_the_child_sample)
    # The goal was to achive at least O(n) complexity, where n is a number of nodes is the tree.
    def edges(self) -> List[Tuple[str, str]]:
        edges = []
        iterate_edges(self.root, final_list=edges)
        return edges
    
    # Methods that returns list containing new trees. Each tree consists of the samples from the selected country.
    # The goal was to achive at least O(n) complexity, where n is a number of nodes is the tree.
    def filter(self, country: str) -> List['Tree']:
        return iterate_country(self.root, country=country)


'''
This function allows recursion to be isolated from the methods in the class (implementation of the Tree.edges method)
It runs in O(n) time and modifies the list passed to it.
'''
def iterate_edges(node: Node, final_list: List):
    if node.children:
        for i in node.children:
            # Adding the node-child edges to the list 
            final_list.append((node.sample.sample_id, i.sample.sample_id))
            #Recursive function calls on children
            iterate_edges(i, final_list)


# This function allows recursion to be isolated from the methods in the class (implementation of the Tree.filter method).
def iterate_country(node: Node, country: str,last_ancestor=None):
    trees = []
    # First checking if the country of the node is the requested country. If not, the ancestor from this function call is being memorised.
    if node.sample.country == country:
        # Else, if it is a suitable country, a new node is created based on a sample of the current node that will belong to the new tree.
        new_node = Node(node.sample)
        # If any ancestor is defined, a new node is added to its children.
        if last_ancestor:
            last_ancestor.children.append(new_node)
        # If not, then the new tree is created and added to the trees list.      
        else:
            trees.append(Tree(new_node))
        # Additionally, the ancestor for recursive function calls is updated.
        last_ancestor_for_future_children = new_node
       
    else: last_ancestor_for_future_children = last_ancestor

    # If a node is not a leaf then the recursive function is called on all of its children
    if node.children:
        for i in node.children:   
            trees += iterate_country(i, country=country, last_ancestor=last_ancestor_for_future_children)
            
    return trees               


# Calculating the edit distance is the function that probably takes the most time. Its pessimistic complexity is theoretically O(m^2).
def edit_disctance(sequence_1, sequence_2):
    # Easy global alignment with mismatch and gap penalties set to -1 to get the edit distance.
    l1 = len(sequence_1)
    l2 = len(sequence_2)
    main_matrix = np.zeros((l1+1,l2+1))
    match_checker_matrix = np.zeros((l1,l2))
    mismatch_penalty = -1
    gap_penalty = -1
    for i in range(l1):
        for j in range(l2):
            if sequence_1[i] != sequence_2[j]: match_checker_matrix[i][j]= mismatch_penalty
    for i in range(l1+1):
        main_matrix[i][0] = i*gap_penalty
    for j in range(l2+1):
        main_matrix[0][j] = j*gap_penalty
    for i in range(1,l1+1):
        for j in range(1,l2+1):
            main_matrix[i][j] = max(main_matrix[i-1][j-1]+match_checker_matrix[i-1][j-1],
                                    main_matrix[i-1][j]+gap_penalty,
                                    main_matrix[i][j-1]+ gap_penalty)
    return -main_matrix[l1, l2]

'''
The function reads sample information from a FASTA file with the given name, creates Sample objects and returns them, sorted in ascending order by sample date.
The goal was to achive O(n * m + n log n) complexity, where n is the number of samples and m is the maximum length of the sequence.
A FASTA file consists of multiple sequences preceded by descriptions:
>sample_ID|country|date
For example: >MT066156.1|Italy|2020-3-13
'''
def read_data(filename: str) -> List[Sample]:
    read = SeqIO.parse(filename, "fasta")
    records = []
    for i in read:
        description = i.description.split('|')
        year, month, day = description[2].split("-")
        sample = Sample(id=description[0], country=description[1],
                        date=datetime.date(year=int(year), month=int(month), day=int(day)), seq=str(i.seq))
        records.append(sample)
    records.sort(key=lambda x: x.date)
    return records

'''
This function returns the optimal phylogenetic tree for the given list of samples. The goal was to achive O(n^2 * m^2) complexity, where n is the number of samples
 and m is the maximum length of the sequence.
The construction of an optimal tree is not complicated and is based on iterating through the samples and checking
the edit distance from a given sample to each older sample and choosing the best one. In this way,
by always choosing the lowest values, we get the optimal tree.
'''
def construct_optimal_tree(samples: List[Sample]):
    samples=[Node(i) for i in samples]
    tree = Tree(root=samples[0])
    for sample in range(1, len(samples)):
        best_match = None
        best_match_score = float("-inf")
        for potential_parent in samples[:sample]:
            # checking if the samples do not have the same date
            if potential_parent.sample.date == samples[sample].sample.date:
                pass
            else:
                ed = edit_disctance(samples[sample].sample.seq, potential_parent.sample.seq)
                if -ed>best_match_score:
                    best_match_score = -ed
                    best_match = potential_parent
        best_match.children.append(samples[sample])
    return tree

'''
A function needed to create an approximate tree, looking for the potential best match among children. 
Note that if there are >1 best matches among the children, the first child found is selected
'''
def best_match_from_chidren(unallocated_node: Node, current_node: Node):
    best_match = None
    best_match_score = float("-inf")
    if current_node.children:
        for child in current_node.children:
            ed = edit_disctance(unallocated_node.sample.seq, child.sample.seq)
            if -ed > best_match_score:
                best_match_score = -ed
                best_match = child
    return best_match_score, best_match


def go_down_the_branch(unallocated_node: Node, branch_highest_node: Node):
    # (a branch is a subtree of the main tree whose root is a child of the root of the main tree)

    dont_stop = True
    potential_parent = branch_highest_node
    best_pp_score = -edit_disctance(unallocated_node.sample.seq, branch_highest_node.sample.seq )
    while dont_stop:
        match_score, match = best_match_from_chidren(unallocated_node=unallocated_node, current_node=potential_parent)
        if match_score >= best_pp_score:
            potential_parent = match
            best_pp_score = match_score
        else:
            dont_stop = False

    return potential_parent, best_pp_score

'''
This function returns the approximate phylogenetic tree for the given list of samples.
The goal was for the function to run noticeably faster than construct_optimal_tree and to return good quality trees.

The algorithm depends on the parameters, the data entered and some random factor.
For each sample, we try to find the sample in the tree with the smallest edit distance from it. When we find it, we add the current sample as a child of the found sample.
We start from the root. The edit distance from it is checked. Then the edit distance from its children is checked. 
The mandatory_check_branches parameter determines what happens when the root has the best edit distance for the sample. 
If it is set to False, the sample is added to the root. If True, the best matching branch + *limit_of_additional_draws* number of random branches will be checked.
Checking a branch is searching through the subtree, looking for the best match. When checking subsequent children, we only go to the one 
that has a better or equal edit distance. If there is no such one, we stop. That is, we simply go down the branch of the tree until we find the best match in it.
If limit_of_additional_draws > 0, then in addition to the best branch, a corresponding number of other randomly selected branches are also checked.
'''
def construct_approximate_tree(samples: List[Sample], mandatory_check_branches=False, limit_of_additional_draws=2):

    samples = [Node(i) for i in samples]
    tree = Tree(root=samples[0])
    root = samples[0]

    for i in range(1, len(samples)):
        sample = samples[i]

        best_matching_branch = None
        best_matching_branch_score = float("-inf")

        branches_visited = []
        # the list above will contain branches that have been visited and are NOT best matches

        # finding the best match among the root children
        for root_child in root.children:
            ed = -edit_disctance(root_child.sample.seq, sample.sample.seq)
            if ed > best_matching_branch_score:
                if best_matching_branch: branches_visited.append(best_matching_branch)
                best_matching_branch = root_child
                best_matching_branch_score = ed
            else:
                branches_visited.append(best_matching_branch)

        # At this point, branches_visited list has all branches except the best one.

        root_score = -edit_disctance(root.sample.seq, sample.sample.seq)
        if root_score > best_matching_branch_score and not mandatory_check_branches:
            root.children.append(sample)
        else:
            if limit_of_additional_draws >= len(branches_visited):
                branches = branches_visited
            else:
                branches = random.sample(branches_visited, limit_of_additional_draws)

            # There might not be any branch in the first iteration
            
            if best_matching_branch: branches.append(best_matching_branch)

            if root_score > best_matching_branch_score:
                best_matching_node = root
                best_matching_node_score = root_score
            else:
                best_matching_node = None
                best_matching_node_score = float("-inf")

            for branch in branches:
                match, score = go_down_the_branch(unallocated_node=sample, branch_highest_node=branch)
                if score > best_matching_node_score:
                    best_matching_node_score = score
                    best_matching_node = match

            best_matching_node.children.append(sample)

    return tree

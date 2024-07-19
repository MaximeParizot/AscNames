from ete3 import Tree
from ete3.coretype.tree import TreeError
import argparse  


#This file compute 3 score to compare ortology relationships prediction models with the 3 groupes tunicates, vertebebrates and outgroup
# First score is the monophyly score corresponding to the group's tendency to be together and without leaf frome other groups 
# Second is the position score corresponding to the relative distance tunicate-vertebrate and tunicate-outgroups , 
# the more the tunicate are closer to vertebrate than the outgroup the higher the score will be 
# Third is the Robison foulds distance between the Tree and a reference Tree who got the same leaf but ina  classic phylogeny, 
# all tunicate together , all vertebrate together and both together far from the outgroup 


#Bad trees are trees with some unresolved leaf so those leaf are juste remove and the score are computed
#Very BAd trees are trees with unresolved nodes so score can't be computed

#Reading arguments
parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="OrthologyRelationships file ",type=str,required=True)
args=parser.parse_args()
file=args.file



def build_reference_tree(groups):
    # Extract groups from dictionary
    equino = groups.get("Equino", [])
    brachio = groups.get("Brachio", [])
    tuniciers = groups.get("tuniciers", [])
    vertebres = groups.get("vertebres", [])

    # Create an empty Tree and root it with the outgroup (Equino)
    reference_tree = Tree(name="root")
    node_equino = reference_tree.add_child(name="equino")

    # Add Equino species to the root node
    for species in equino:
        node_equino.add_child(name=species)

    # Create a node for brachio and tuniciers_vertebres
    brachio_tuniciers_vertebres = reference_tree.add_child(name="brachio_tuniciers_vertebres")

    # Add a node for brachio
    node_brachio = brachio_tuniciers_vertebres.add_child(name="brachio")
    for species in brachio:
        node_brachio.add_child(name=species)

    # Create a subtree for tuniciers and vertebres
    tuniciers_vertebres = brachio_tuniciers_vertebres.add_child(name="tuniciers_vertebres")

    # Add a node for tuniciers
    node_tuniciers = tuniciers_vertebres.add_child(name="tuniciers")
    for species in tuniciers:
        node_tuniciers.add_child(name=species)

    # Add a node for vertebres
    node_vertebres = tuniciers_vertebres.add_child(name="vertebres")
    for species in vertebres:
        node_vertebres.add_child(name=species)

    return reference_tree

#Compute the monophyly score for a group 
def monophyly_proportion(tree, group):
    #print(tree)
    if len(group)==1: 
        return 1
    #Find the Most Recent Commun Ancestor (MRCA) of the group in the tree
    mrca = tree.get_common_ancestor(group)
    group_leaves = set(mrca.get_leaf_names())
    #Matched leaves contains the leaf that are above the mrca and in the group 
    matched_leaves = group_leaves.intersection(set(group))
    #the monophyly score is the proportions of leaf above the mrca that are in the group
    return len(matched_leaves) / len(group_leaves)

#Compute the position score of a Tree
def relative_position(tree, outgroup, group1, group2):
    #group1 : tunicates ; group 2 : vertebrates
    #find mrca of for each groups 
    outgroup_mrca = tree.get_common_ancestor(outgroup)
    group1_mrca = tree.get_common_ancestor(group1)
    group2_mrca = tree.get_common_ancestor(group2)

    #Compute the distance bewteeen the mrca of each group 
    dist_group1_group2 = group1_mrca.get_distance(group2_mrca)
    dist_group1_outgroup = group1_mrca.get_distance(outgroup_mrca)
    
    if outgroup_mrca==tree.get_tree_root() : 
        dist_group1_outgroup=2*dist_group1_outgroup
    # Verify if all the distance are null to avoid division by zero
    if dist_group1_group2 == 0 and dist_group1_outgroup == 0 :
        return 0.0  # Retourner 0 si l'une des distances est nulle
    
    # Compute the position score 
    relative_pos = 1-((dist_group1_group2) / 
                    (dist_group1_group2 + dist_group1_outgroup))
    
    
    return relative_pos

#Compute the Robison Foulds distance between two trees
def calculate_rf_distance(tree1, tree2):
        results = tree1.robinson_foulds(tree2)
        #results contain a lot of things we just want the rf and the max rf
        rf, max_rf = results[0], results[1]
        #then normalized it to get a distance 
        normalized_rf = rf / max_rf
        return normalized_rf

#Compute Composite score 
def calculate_composite_score(observed_tree, reference_tree, groups):
    monophyly_score = 0
    sgroups={'vertebres': groups['vertebres'],'tuniciers' : groups['tuniciers']}
    max_score = len(sgroups)  # Number of groups 

    # compute monophyly score for each group in the observed tree
    for group_name, group_species in sgroups.items():
        monophyly_score += monophyly_proportion(observed_tree, group_species)
    
    # then the position score of the tree 
    position_score=relative_position(observed_tree, groups["Equino"], groups["tuniciers"], groups["vertebres"])
    # then the rf distance
    rf_distance = calculate_rf_distance(observed_tree, reference_tree)

    # Normalized the monphyly score 
    monophyly_score = monophyly_score / max_score

    # Then compute the composite score which is the mean of those 3 score 
    composite_score = (monophyly_score + position_score + (1 - rf_distance)) / 3

    return monophyly_score,position_score,1-rf_distance,composite_score
    #print(f"monophyly_score: {monophyly_score:.2f}")
    #print(f"position_score: {position_score:.2f}")
    #print(f"rf_distance): {rf_distance:.2f}")
    #return composite_score


def remove_unresolved_nodes(tree):
    root = tree.get_tree_root()
    unresolved_nodes = []
    
    # Look at each child of the root 
    for child in root.get_children():
        if child.is_leaf() :
            unresolved_nodes.append(child)
    
    # Remove unresolved node if there are more than 2 children 
    for node in unresolved_nodes:
        if len(root.get_children()) > 2:
            node.detach()
    
    return tree

# Verify if a tree is rooted (= root get 2 and only 2 children)
def is_rooted(tree):
    root = tree.get_tree_root()
    return len(root.get_children()) == 2



def contains_all_groups(node, tuniciers, vertebres, equino):
    # Check if a node or its descendants contain species from all three groups
    species = node.get_leaf_names()
    has_tuniciers = any(species_name.split('|')[0] in tuniciers for species_name in species)
    has_vertebres = any(species_name.split('|')[0] in vertebres for species_name in species)
    has_equino = any(species_name.split('|')[0] in equino for species_name in species)
    return has_tuniciers and has_vertebres and has_equino


def process_tree(node, tuniciers, vertebres, equino, scores):
    #Apply scoring method to the good trees 
    subtree_contains_all = False

    #For all child of current node wich coutains one tunicate one vertebrate and one equino, continu the processus with those children
    for child in node.children:
        if contains_all_groups(child, tuniciers, vertebres, equino):
            process_tree(child, tuniciers, vertebres, equino, scores)
            subtree_contains_all = True

    #If there are none score the tree rooted in the current node
    if not subtree_contains_all:
        Brachio=['Bbel','Blan']
        leaf_names = node.get_leaf_names()
        #Define groups 
        groups = {
            "Equino": [name for name in leaf_names if name.split('|')[0] in equino],
            "tuniciers": [name for name in leaf_names if name.split('|')[0] in tuniciers],
            "vertebres": [name for name in leaf_names if name.split('|')[0] in vertebres],
            "Brachio": [name for name in leaf_names if name.split('|')[0] in Brachio]
        }
        #build the reference tree
        reference_tree=build_reference_tree(groups)
        #Compute scores 
        subtree_score = calculate_composite_score(node,reference_tree,groups)
        for i in range(0,len(subtree_score)) : 
            scores['scores'][i]+=subtree_score[i]
        scores['nbtrees'] += 1

def analyze_tree(tree):
    #Initialize processus of cuting/scoring to a tree
    #Define groups 
    equino = ['Spur','Apla','Pmin']
    tuniciers =['Phmamm','Phfumi','Cisavi','Cirobu','Moocci','Moocul','Mooccu','Boschl','Boleac','Haaura','Harore','Coinfl','Stclav']
    vertebres = ['Cmil','Lcha','Hsap','Mmus','Ggal','Psin']
    #Initialize scores and numbers a subtree scored in the tree , then return them 
    scores = {'scores': [0,0,0,0], 'nbtrees': 0}
    process_tree(tree, tuniciers, vertebres, equino, scores)
    return scores['scores'], scores['nbtrees']


mnscoremoyen=0
pstscoremoyen=0
rfdistmoyen=0
cpscoremoyen=0
badtrees=0
verybadtrees=0
nm=0

#Score all the Tree in a file 
with open(file) as f: 
    for line in f.readlines() : 
        print(nm)
        observed_tree = Tree(line)
        #Try to score the tree, if an error occure it's because it's a bad tree (coutains unresolved nodes for example)
        try : 
            scores,nbtree=analyze_tree(observed_tree)
            nm+=nbtree
            mnscoremoyen+=scores[0]
            pstscoremoyen+=scores[1]
            rfdistmoyen+=scores[2]
            cpscoremoyen+=scores[3]
        except TreeError :
            badtrees+=1
    
#Print the mean of each score and number of bad trees
print("Score moyens : ")
print(f"Composite Score: {cpscoremoyen/nm:.4f}")
print(f"Monophyly Score: {mnscoremoyen/nm:.4f}")
print(f"Position Score: {pstscoremoyen/nm:.4f}")
print(f"RF score: {rfdistmoyen/nm:.4f}")
print("Nb bad trees : "+str(badtrees)+'/'+str(nm+badtrees))
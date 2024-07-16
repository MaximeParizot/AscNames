from ete3 import Tree
from ete3.coretype.tree import TreeError
import argparse  

#Reading arguments
parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="OrthologyRelationships file ",type=str,required=True)
args=parser.parse_args()
file=args.file


def build_reference_tree(groups):
    # Extract groups from dictionnary
    outgroup = groups.get("outgroup", [])
    tuniciers = groups.get("tuniciers", [])
    vertebres = groups.get("vertebres", [])

    # Create an empty Tree
    reference_tree = Tree()

    # Add each group nodes 
    node_outgroup = reference_tree.add_child(name="outgroup")
    node_tuniciers = reference_tree.add_child(name="tuniciers")
    node_vertebres = reference_tree.add_child(name="vertebres")

    # Add each species in the corresponding groups 
    for species in outgroup:
        node_outgroup.add_child(name=species)

    for species in tuniciers:
        node_tuniciers.add_child(name=species)

    for species in vertebres:
        node_vertebres.add_child(name=species)

    # Root the tree with the outgroup
    reference_tree.set_outgroup(node_outgroup)

    return reference_tree

#Compute the monophyly score for a group 
def monophyly_proportion(tree, group):
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
    dist_group2_outgroup = group2_mrca.get_distance(outgroup_mrca)
    
    # Verify if all the distance are null to avoid division by zero
    if dist_group1_group2 == 0 and dist_group1_outgroup == 0 and dist_group2_outgroup == 0:
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
    max_score = len(groups)  # Number of groups 

    # compute monophyly score for each group in the observed tree
    for group_name, group_species in groups.items():
        monophyly_score += monophyly_proportion(observed_tree, group_species)

    # then the position score of the tree 
    position_score=relative_position(observed_tree, groups["outgroup"], groups["tuniciers"], groups["vertebres"])
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



# Define groups 
outgroup = ['Spur','Apla','Bbel','Blan','Ajap','Pmin']
tuniciers =['Phmamm','Phfumi','Cisavi','Cirobu','Moocci','Moocul','Mooccu','Boschl','Boleac','Haaura','Harore','Coinfl','Stclav']
vertebres = ['Cmil','Lcha','Hsap','Mmus','Ggal','Psin']

mnscoremoyen=0
pstscoremoyen=0
rfdistmoyen=0
cpscoremoyen=0
badtrees=0
verybadtrees=0

# Look at each line(=Tree) of the file 
with open(file) as f: 
    nm=0
    for line in f.readlines():
        nm+=1
        print(nm)
        observed_tree = Tree(line)
        #Check if the tree is rooted and if not try to remove the unresolved nodes 
        if not(is_rooted(observed_tree)): 
            observed_tree=remove_unresolved_nodes(observed_tree)
            badtrees+=1
        leaf_names = observed_tree.get_leaf_names()
        #Define the group of the observed and reference tree 
        groups = {
            "outgroup": [name for name in leaf_names if name.split('|')[0] in outgroup],
            "tuniciers": [name for name in leaf_names if name.split('|')[0] in tuniciers],
            "vertebres": [name for name in leaf_names if name.split('|')[0] in vertebres]
        }

        #build the reference tree
        reference_tree = build_reference_tree(groups)
        #then compute all the score 
        try : 
            monophyly_score,position_score,rf_distance,composite_score = calculate_composite_score(observed_tree, reference_tree, groups)
            mnscoremoyen+= monophyly_score
            pstscoremoyen+=position_score
            rfdistmoyen+=rf_distance
            cpscoremoyen+=composite_score
        #If tree is still unresolved, don't look at it 
        except TreeError : 
            #print(tree)
            verybadtrees+=1


print("Mean scores  : ")
print(f"Composite Score: {cpscoremoyen/nm:.4f}")
print(f"Monophyly Score: {mnscoremoyen/nm:.4f}")
print(f"Position Score: {pstscoremoyen/nm:.4f}")
print(f"RF score: {rfdistmoyen/nm:.4f}")
print("Nb bad trees : "+str(badtrees)+'/'+str(nm))
print("Nb very bad trees : "+str(verybadtrees)+'/'+str(nm))
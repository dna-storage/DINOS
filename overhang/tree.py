'''
Author: Kevin Volkel

Filename: tree.py

Description: Provides tools that analyze a set of strands and generates a directed acyclic graph from the set of strands

'''
#import networkx as nx #graph support library
import math
import sys
from overhang.util.overhang_utils import * #get/cut overhang utility functions here
import overhang.dnastorage_utils.codec.base_conversion as bc
from reaction_node import *

#Steps to correcting a node
'''
    1. check to see if the node can be merged with its grand child,
       - node can be merged if it has 1 real child, and that child has 1 real child, e.g. 1 real grandchild
       - the node must be a grandchild since that means they both are using the same overhangs to merge substrands
       - when merging, remove the child node, and all of its padding children, remove the initial node and move the grandchild along with all of its padding up to the original node's spot
    2. If the node was merged, leave, the node moved up should already have its padding

    3. If the node was not merged, padding should be added 
       - padding is added based on the different strategies employed
        Strategy 1 (ID'd as 0): We have m^2 NOP strands, which means a NOP exists for every overhang permutation, meaning only one NOP shall be inserted to fix the strand's overhangs. To figure out the NOP strand, find the beginning overhang of the input strand get the ID, get the end overhang and get its ID. The end ID should be (beg_ID+(m-1))%(m-1) [if indexing starting from 0] or (beg_ID+1)%(m-1) : if the height is either odd or even respectively (assuming height counts start at 1). After Indentifying the end ID, the NOP strand can be derived as old_end:NOP:new_end. Add this strand to the node's padding field.
      
       Strategy 2 (ID'd as 1): We have 2m NOP strands, which means we have m NOPS for both cases, a even height and an odd height. Same as strategy 1, ID the start overhang, the old end and the new end. calculate how far the old end is from the new end to figure out how many NOP strands to insert into the node, which should be between 1 and (m-1)-1. Calculate the concantenation of all of the NOPs added, give this to the input node and return

       Strategy 3 (ID'd as 2): We have m NOP strands, which means we have m Nops only for the case of the odd height. To derive appropriate strands for the even height, we need to add in new child nodes to the input nodes. Steps are same as strategy 2, where we find out the NOP strands to insert as if we had the 2m strands, but if the height is even, we must convert these NOPs to new reactions who join together (m-1) NOP strands from the original set of NOP strands. The padding strand given to this node should reflect the (m-1)*R NOP strands that need to be inserted. The reactions that create these padding strands are shareable.

   4. After adding the padding, and counting the padding nodes, the node should be updated and the method can be returned from.


NOTE: The ultimate strand output from a node is going the be the concantenation of each of its children outputs, along with the concantenation of the padding strand. If no padding is needed anywhere in the tree, then the concantenation (making sure not to double count overhangs) should be the same as the originally calculated output for the node. 

'''
def node_correction(tree,node,parent):
    if node.get_true_children()==1 and node.get_child_list()[0].get_true_children()==1:
        #print "merging nodes start at height {}, total node count {} iteration {}".format(node.get_height(),tree.order(),node_correction.counter)
        #print node.get_child_list()[0].get_child_list()[0].get_height()
        # we can merge the current node with its grand child
        remove_count=1+len(node.get_child_list()[0].get_child_list()) #remove the current node, the child node, and all of its children (not including the child that is the one we merging)
        #print remove_count
        tree.dec_node_count(remove_count) #record removal in the tree
        tree.dec_pad_count(len(node.get_child_list())-node.get_true_children()) #decrement non-true children 
        #since we traverse depth first, last entry in the parent list should be the current node
        parent.get_child_list()[-1]=node.get_child_list()[0].get_child_list()[0] #grandchild should also already be padded
    else:#can't merge, so we need to pad the node out
        strand=""
        strand=node.get_corrected_substrand()
        #print node.get_height()
        #print strand
        
        start_overhang=get_start_overhang(strand,tree.overhang_length)
        end_overhang=get_end_overhang(strand,tree.overhang_length)
        startID=tree.overhang_to_num[start_overhang]
        endID=tree.overhang_to_num[end_overhang]
        if node.get_height()&0x000000001:
            #node is odd height: final end ID should be (startID+(m-1))%(m)
            final_endID=(startID+tree.num_overhangs-1)%(tree.num_overhangs)
            if(final_endID>endID):
                end_difference=final_endID-endID
            elif(final_endID<endID): #wrap around situation
                distance_to_wrap=(tree.num_overhangs-1)-endID #find how far from left hand of wrap
                distance_from_wrap=final_endID+1 #find how far from right hand of wrap
                end_difference=distance_from_wrap+distance_to_wrap
            assert (end_difference>0 and end_difference<tree.num_overhangs-1)
        else:
            #print "fixing even node"
            #print "end {} {}".format(endID,end_overhang)
            #print "start {} {}".format(startID,start_overhang)
            #node is even height: final end ID should be (startID+1)%(m+1)
            final_endID=(startID+1)%(tree.num_overhangs)
            assert(not startID==endID)
            if(final_endID<endID):
                end_difference=endID-final_endID
            elif (final_endID>endID):
                distance_to_wrap=endID+1
                distance_from_wrap=(tree.num_overhangs-1)-final_endID
                end_difference=distance_from_wrap+distance_to_wrap
            assert (end_difference>0 and end_difference<tree.num_overhangs-1)

        final_end_overhang=tree.num_to_overhang[final_endID]
        #at this point, we know what the end overhang should be, what it is, and how many NOPs should be inserted if we are using an expanded direct repair technique (1)
        
        if tree.end_repair_technique==0:
            #m^2 repair: direct repair
            repair_strand=end_overhang+'NULL'+final_end_overhang
            node.set_pad_strand(repair_strand) #repair strand will just be the current end plus the final end
            return 
        elif tree.end_repair_technique==1:
            #2m repair: expanded direct repair
            if node.get_height()&0x00000001:
                #node is odd height: increment from the current endID
                repair_strand=end_overhang
                new_endID=endID
                while(new_endID!=final_endID):
                    repair_strand+='NULL'
                    new_endID=(new_endID+1)%(tree.num_overhangs)
                    repair_strand+=tree.num_to_overhang[new_endID]
            else:
                #node is even height: decrement from current endID
                repair_strand=end_overhang
                new_endID=endID
                while(new_endID!=final_endID):#create pad strand while we haven't reached the right end
                    repair_strand+='NULL'
                    new_endID=(new_endID-1)
                    if(new_endID<0): new_endID=tree.num_overhangs-1
                    repair_strand+=tree.num_to_overhang[new_endID]
            node.set_pad_strand(repair_strand)
            return
        elif tree.end_repair_technique==2:
            #m repair: additional node repair
            #nodes only need to be added if this node is even
            if node.get_height()&0x00000001:
                #node is odd, treat it like technique 1
                repair_strand=end_overhang
                new_endID=endID
                while(new_endID!=final_endID):
                    repair_strand+='NULL'
                    new_endID=(new_endID+1)%(tree.num_overhangs)
                    repair_strand+=tree.num_to_overhang[new_endID]
                node.set_pad_strand(repair_strand)
                return
            else: 
                #node is even, need to generate a reaction for each "technique 1" strand since we do not have them directly already
                new_endID=endID
                repair_strand=""
                while(new_endID!=final_endID):
                    reaction_start=new_endID
                    new_endID=(new_endID-1)
                    if(new_endID<0): new_endID=tree.num_overhangs-1
                    reaction_end=new_endID
                    #create the strand that would output from the padding reaction
                    reaction_strand=tree.num_to_overhang[reaction_start]
                    reaction_iter=0
                    while(reaction_start!=reaction_end):
                        reaction_strand+='NULL'
                        reaction_start=(reaction_start+1)%(tree.num_overhangs)
                        reaction_strand+=tree.num_to_overhang[reaction_start]
                        reaction_iter+=1
                    assert(reaction_iter==tree.num_overhangs-1) #loop should complete a full reaction
                    if reaction_strand in tree.strand_hash_table:
                        #print "reuse padding"
                        padNode=tree.strand_hash_table[reaction_strand]
                        assert(padNode.is_pad())
                        padNode.inc_use_count()
                        padNode.insert_strand_ID(node._strandIDList)#inherit the node's ID list
                    else:
                        #print "new node"
                        padNode=ReactionNode(node._strandIDList,-1,True) #make pad reaction
                        tree.inc_node_count()
                        tree.inc_pad_count()
                        tree.strand_hash_table[reaction_strand]=padNode
                        
                    node.add_child(padNode)#add padding node to the current node
                    assert(len(node.get_child_list())<tree.num_overhangs)
                    repair_strand+=cut_end_overhang(reaction_strand,tree.overhang_length) #cut the end overhang off, accumulating the final repair strand here
                repair_strand+=tree.num_to_overhang[new_endID]#add final overhang
                node.set_pad_strand(repair_strand)
                return






def ideal_tree_rec(s,parent_node,h,strand_index,reactiontree):
    block_group_size=(reactiontree.num_overhangs-1)**h
    for i in range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)):#breaks up the input strand into sub reactions  
        if i == (len(s)-(reactiontree.overhang_length)): break
        end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
        #assert end_index<=len(s)
        if end_index>=len(s): end_index=len(s)
        substrand=s[i:end_index]
        datastrand=""
        #cut data out of substrand
        for data_start in range(reactiontree.overhang_length,len(substrand),reactiontree.codeword_length+reactiontree.overhang_length):
            datastrand+=substrand[data_start:data_start+reactiontree.codeword_length]
        if datastrand in reactiontree.strand_hash_table and reactiontree.strand_hash_table[datastrand].get_height()==h:
            if reactiontree.h_array is not None:
                reactiontree.h_array[h-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
            node=reactiontree.strand_hash_table[datastrand]
            node.inc_use_count()
            node.insert_strand_ID(strand_index)
            if parent_node!=None:
                #reactiontree.add_edge(node,parent_node)
                parent_node.add_child(node)
            continue #move onto next substrand at the same level without traversing deeper down the tree
        else:
            node=ReactionNode(strand_index,h,False)
            reactiontree.strand_hash_table[datastrand]=node #add node to the hastable
            reactiontree.inc_node_count()
            if parent_node!=None:
                #reactiontree.add_edge(node,parent_node)
                parent_node.add_child(node)
            else: #this is a terminator
                reactiontree.add_terminator(node)
                node.set_term()
            #before going to the next sub-strand at height h, move down the tree
            if h >1: #if we are at height 1, we are at the bottom of the tree
                ideal_tree_rec(substrand,node,h-1,strand_index,reactiontree)
                node.set_true_children(len(node.get_child_list()))
    return


#optimal tree construction that takes advantage of any data repeat, but also constructs the tree recursively as to not double count sub reactions of a matching strand
def construct_tree_ideal(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    for strand_index, s in enumerate(strands):
        h=int(math.log(strand_length_in_codewords,num_overhangs-1))
        ideal_tree_rec(s,None,h,strand_index,reactiontree)
    return reactiontree




#recursive method used to build the base-opt tree
def baseopt_tree_rec(s,parent_node,h,strand_index,reactiontree):
    block_group_size=(reactiontree.num_overhangs-1)**h
    for i in range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)):#breaks up the input strand into sub reactions  
        if i == (len(s)-(reactiontree.overhang_length)): break
        end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
        if end_index >= len(s): end_index=len(s)
        #assert end_index<=len(s)
        substrand=s[i:end_index]
        if substrand in reactiontree.strand_hash_table and reactiontree.strand_hash_table[substrand].get_height()==h: #make sure heights match
            if reactiontree.h_array is not None:
                reactiontree.h_array[h-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
            node=reactiontree.strand_hash_table[substrand]
            node.inc_use_count()
            node.insert_strand_ID(strand_index)
            if parent_node!=None:
                #reactiontree.add_edge(node,parent_node)
                parent_node.add_child(node)
            continue #move onto next substrand at the same level without traversing deeper down the tree
        else:
            node=ReactionNode(strand_index,h,False)
            reactiontree.strand_hash_table[substrand]=node #add node to the hastable
            reactiontree.inc_node_count()
            if parent_node!=None:
                #reactiontree.add_edge(node,parent_node)
                parent_node.add_child(node)
            #before going to the next sub-strand at height h, move down the tree
            else:
                #this is a terminator
                reactiontree.add_terminator(node)
                node.set_term()
            if h >1: #if we are at height 1, we are at the bottom of the tree
                baseopt_tree_rec(substrand,node,h-1,strand_index,reactiontree)
                node.set_true_children(len(node.get_child_list()))#should be done with all sub tree work once reaching here
                c_strand=create_corrected_strand_short(node,substrand,reactiontree.overhang_length)
                node.set_corrected_substrand(c_strand)
                if node.get_true_children()<reactiontree.ideal_children_nodes and not node.is_term() and node.get_pad_strand() is "":
                    node_correction(reactiontree,node,parent_node)
                elif node.is_term():
                    node.set_corrected_substrand("")#remove term's strand
                for _ in node.get_child_list():
                    _.set_corrected_substrand("") #no longer need these children's corrected substrands
                
            if h==1: #at the base but we still need to check for node correction
                node.set_corrected_substrand(substrand) #bottom of tree set corrected to current
                start_overhang=get_start_overhang(substrand,reactiontree.overhang_length)
                #print "node strand {}".format(substrand)
                #print "start overhang {}".format(start_overhang)
                end_overhang=get_end_overhang(substrand,reactiontree.overhang_length)
                #print "end overhang {}".format(end_overhang)
                start_overhangID=reactiontree.overhang_to_num[start_overhang]
                end_overhangID=reactiontree.overhang_to_num[end_overhang]
                #print"start ID {} end ID {}".format(start_overhangID,end_overhangID)
                #end ID should be m-1 larger (mod m) than the start overhang
                ideal_ID=(start_overhangID+(reactiontree.num_overhangs-1))%(reactiontree.num_overhangs)
                if not ideal_ID==end_overhangID: #end repair needs to be done at the base
                    node_correction(reactiontree,node,parent_node)
                node.set_true_children(reactiontree.num_overhangs-1)
                
    return


def construct_tree_baseopt(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    for strand_index, s in enumerate(strands):
        h=int(math.log(strand_length_in_codewords,num_overhangs-1))
        baseopt_tree_rec(s,None,h,strand_index,reactiontree)
    return reactiontree
        




#lite version of ideal recursive tree build
def ideal_tree_rec_lite(s,h,strand_index,reactiontree):
    block_group_size=(reactiontree.num_overhangs-1)**h
    for i in range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)):#breaks up the input strand into sub reactions
        if i == (len(s)-(reactiontree.overhang_length)): break
        end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
        #assert end_index<=len(s)
        if end_index>=len(s): end_index=len(s)
        substrand=s[i:end_index]
        datastrand=""
        #cut data out of substrand
        for data_start in range(reactiontree.overhang_length,len(substrand),reactiontree.codeword_length+reactiontree.overhang_length):
            datastrand+=substrand[data_start:data_start+reactiontree.codeword_length]
        #if h==12:
            #print datastrand
            #print substrand
            #print len(datastrand)
        if datastrand in reactiontree.strand_hash_table and reactiontree.strand_hash_table[datastrand]>>32==h:
            if reactiontree.h_array is not None:
                reactiontree.h_array[h-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
            use_count=reactiontree.strand_hash_table[datastrand]&(4294967296-1)
            use_count+=1
            height=reactiontree.strand_hash_table[datastrand]>>32
            reactiontree.strand_hash_table[datastrand]=(height<<32)|use_count
            continue #move onto next substrand at the same level without traversing deeper down the tree
        else:
            hash_table_value=(h<<32)|1 #height in top 32, use count in bottom 32
            reactiontree.strand_hash_table[datastrand]=hash_table_value
            reactiontree.inc_node_count()
            #before going to the next sub-strand at height h, move down the tree
            if h >1: #if we are at height 1, we are at the bottom of the tree
                ideal_tree_rec_lite(substrand,h-1,strand_index,reactiontree)
    return


#entry point to extremely lightweight node counting for the base optimzation case, also takes into account padding needed for incomplete trees
def construct_tree_ideal_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    for strand_index, s in enumerate(strands):
        h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
        ideal_tree_rec_lite(s,h,strand_index,reactiontree)
    #create a unoptimized tree lite to calculate necessary padding overheads
    unopt_tree=construct_tree_unoptimized_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy,h_array)
    reactiontree.add_pad_count(unopt_tree.get_pad_count()/len(strands)) #can optimize therefore the price of padding count should be only paid once, padding is also not data dependent only alignment dependent per strand
    reactiontree.add_node_count(unopt_tree.get_pad_count()/len(strands))
    reactiontree.dec_node_count(unopt_tree.get_singleton_nodes_removed())
    return reactiontree






#recursive method, with modifications to make optimization reaction countring more lightweight
def baseopt_tree_rec_lite(s,h,strand_index,reactiontree):
    block_group_size=(reactiontree.num_overhangs-1)**h
    for i in range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)):#breaks up the input strand into sub reactions  
        if i == (len(s)-(reactiontree.overhang_length)): break
        end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
        if end_index >= len(s): end_index=len(s)
        #assert end_index<=len(s)
        substrand=s[i:end_index]
        if substrand in reactiontree.strand_hash_table and reactiontree.strand_hash_table[substrand]>>32==h: #make sure heights match
            if reactiontree.h_array is not None:
                reactiontree.h_array[h-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
            use_count=reactiontree.strand_hash_table[substrand]&(4294967296-1)
            use_count+=1
            height=reactiontree.strand_hash_table[substrand]>>32
            reactiontree.strand_hash_table[substrand]=(height<<32)|use_count
            continue #move onto next substrand at the same level without traversing deeper down the tree
        else:
            hash_table_value=(h<<32)|1 #height in top 32, use count in bottom 32
            reactiontree.strand_hash_table[substrand]=hash_table_value
            reactiontree.inc_node_count()
            if h >1: #if we are at height 1, we are at the bottom of the tree
                baseopt_tree_rec_lite(substrand,h-1,strand_index,reactiontree)
                
    return

#entry point to extremely lightweight node counting for the base optimzation case, also takes into account padding needed for incomplete trees
def construct_tree_baseopt_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    for strand_index, s in enumerate(strands):
        h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
        baseopt_tree_rec_lite(s,h,strand_index,reactiontree)
    #create a unoptimized tree lite to calculate necessary padding overheads
    unopt_tree=construct_tree_unoptimized_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy,h_array)
    reactiontree.add_pad_count(unopt_tree.get_pad_count()/len(strands)) #can optimize therefore the price of padding count should be only paid once, padding is also not data dependent only alignment dependent per strand
    reactiontree.dec_node_count(unopt_tree.get_singleton_nodes_removed()) #account for removal of nodes
    return reactiontree





#recursive method used to build the unoptimized tree
def unopt_tree_rec(s,parent_node,h,strand_index,reactiontree):
    block_group_size=(reactiontree.num_overhangs-1)**h
    for i in range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)):#breaks up the input strand into sub reactions  
        if i == (len(s)-(reactiontree.overhang_length)): break
        end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
        if end_index>=len(s): end_index=len(s)
        #assert end_index<=len(s)
        substrand=s[i:end_index]
        node=ReactionNode(strand_index,h,False)
        reactiontree.inc_node_count()
        if parent_node!=None:
            #reactiontree.add_edge(node,parent_node)
            parent_node.add_child(node)
        #before going to the next sub-strand at height h, move down the tree
        else:
            #this is a terminator
            reactiontree.add_terminator(node)
            node.set_term()
        if h >1: #above the base of the tree
            unopt_tree_rec(substrand,node,h-1,strand_index,reactiontree)
            node.set_true_children(len(node.get_child_list()))#should be done with all sub tree work once reaching here
            c_strand=create_corrected_strand_short(node,reactiontree.overhang_length)
            #print node.get_height()
            node.set_corrected_substrand(c_strand)
            #print node.get_corrected_substrand()
            #print node.get_substrand()
            if node.get_true_children()<reactiontree.ideal_children_nodes and not node.is_term() and node.get_pad_strand() is "":
                node_correction(reactiontree,node,parent_node)
            elif node.is_term():
                node.set_corrected_substrand("")
            for _ in node.get_child_list():
                _.set_corrected_substrand("") #no longer need these children's corrected substrands

        if h==1: #at the base but we still need to check for node correction
            start_overhang=get_start_overhang(substrand,reactiontree.overhang_length)
            end_overhang=get_end_overhang(substrand,reactiontree.overhang_length)
            start_overhangID=reactiontree.overhang_to_num[start_overhang]
            end_overhangID=reactiontree.overhang_to_num[end_overhang]
            #end ID should be m-1 larger (mod m) than the start overhang
            ideal_ID=(start_overhangID+(reactiontree.num_overhangs-1))%(reactiontree.num_overhangs)
            if not ideal_ID==end_overhangID: #end repair needs to be done at the base
                node_correction(reactiontree,node,parent_node)
            node.set_true_children(reactiontree.num_overhangs-1)
            node.set_corrected_substrand(substrand) #bottom of tree set corrected to current
    return

def construct_tree_unoptimized(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    #reactiontree_x=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    '''
    TODO: Based on strand length and the number of overhangs, do quick analytical calculations rather than recursive algorithm
    '''
    reactiontree=analytical_ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy)
    
    #print "analytical model"
    h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
    for strand_index, s in enumerate(strands):
        unopt_tree_rec(s,None,h,strand_index,reactiontree_x)
    
    return reactiontree

def construct_tree_unoptimized_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    #reactiontree_x=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    '''
    TODO: Based on strand length and the number of overhangs, do quick analytical calculations rather than recursive algorithm    '''
    reactiontree=analytical_ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy)
    #print "analytical model"
    h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
    pad_nodes=0
    #print h
    total_nodes=0
    total_inserts=0
    previous=False #flag to indicate the previous node is an only child
    previous_padding=0
    previous_inserts=0
    singleton_nodes_removed=0
    for i in range(1,h+1)[::-1]:
        #go from the tree top down
        if i==h:
            total_nodes+=1
        else:
            width=int(math.ceil(float(strand_length_in_codewords)/float((num_overhangs-1)**i)))
            if i==(h-1):
                total_nodes+=width
                continue
            else:
                if width%(num_overhangs-1)==0:
                    total_nodes+=width
                    previous=0
                else:
                    #need to handle incomplete portion of tree
                    
                    if previous and width%(num_overhangs-1)==1:
                        #print "no-opt removal height {} starting reaction count {}".format(i+2,total_nodes+width)
                        #previous is set,therefore we can remove the previous 2 nodes
                        total_nodes+=(width-1) #add width-1 (don't need last node)
                        total_nodes-=1 #(remove parent node)
                        total_nodes-=previous_padding #(no longer need padding of last node)
                        #print "no-opt nodes after removal {}".format(total_nodes)
                        pad_nodes-=previous_padding
                        total_inserts-=previous_inserts
                        previous=False
                        previous_padding=0
                        singleton_nodes_removed+=2
                    else:
                        total_nodes+=width
                        _pad=0
                        _inserts=0
                        #need to add in some padding
                        if repair_strategy==0:
                            #m^2 direct repair, only need one insert
                            _inserts=1
                        elif repair_strategy==1:
                            #2m direct repair
                            _inserts=(num_overhangs-1)-width%(num_overhangs-1)
                        elif repair_strategy==2:
                            _pad=(num_overhangs-1)-width%(num_overhangs-1)
                            _inserts=_pad*(num_overhangs-1)#each pad node should have m-1 inserts
                            total_nodes+=_pad
                        total_inserts+=_inserts
                        pad_nodes+=_pad
                        if width%(num_overhangs-1)==1:
                            #print "hang off {}".format(i)
                            #indicate a lonely end node in the width
                            previous=True
                            previous_padding=_pad
                            previous_inserts=_inserts
                        else:
                            previous=False
                            previous_padding=0
                            previous_inserts=0
    total_nodes*=len(strands)
    pad_nodes*=len(strands)
    total_singleton_nodes_removed=len(strands)*singleton_nodes_removed
    reactiontree.add_node_count(total_nodes)
    reactiontree.add_pad_count(pad_nodes)
    reactiontree.set_insert_per_strand(total_inserts)
    reactiontree.set_singleton_nodes_removed(total_singleton_nodes_removed)

    return reactiontree






########################## Breadth First Tree Construction (Going to be useful for Assembly Tree Changes for Optimization)#######
#entry point to extremely lightweight node counting for the base optimzation case, also takes into account padding needed for incomplete trees
def construct_tree_baseopt_lite_w(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    for strand_index, s in enumerate(strands):
        #print len(s)
        h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
        skip_array=[False]*strand_length_in_codewords #If an entry is True: you skip this sub strand
        for height in range(1,h+1)[::-1]:
            block_group_size=(reactiontree.num_overhangs-1)**height
            for strandID, i in enumerate(range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length))):#breaks up the input strand into sub reactions
                #print strandID
                strandIDmod=strandID%num_overhangs
                if i == (len(s)-(reactiontree.overhang_length)): break
                if skip_array[strandID*block_group_size]: continue #this tree is already optimized out
                end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
                if end_index>=len(s): end_index=len(s)
                substrand=s[i:end_index]
                datastrand=""
                #cut data out of substrand
                for data_start in range(reactiontree.overhang_length,len(substrand),reactiontree.codeword_length+reactiontree.overhang_length):
                    datastrand+=substrand[data_start:data_start+reactiontree.codeword_length]
                #print len(datastrand)
                #print codeword_length*block_group_size
                assert len(datastrand)==codeword_length*block_group_size
                if datastrand in reactiontree.strand_hash_table:
                    if strandIDmod in reactiontree.strand_hash_table[datastrand]:#check overhang requirement
                        #have a match, take the height statistic, don't count the reaction, mark the skip_array
                        if reactiontree.h_array is not None:
                            reactiontree.h_array[height-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
                        for _ in range(0,block_group_size):
                            skip_array[_+strandID*block_group_size]=True #set the underlying codewords to all True
                        reactiontree.strand_hash_table[datastrand][strandIDmod]+=1
                        continue #move onto next substrand
                    else: #have seen the data, but not this strandIDmod (overhang version), increment the reaction count
                        reactiontree.strand_hash_table[datastrand][strandIDmod]=1
                        reactiontree.inc_node_count()
                        continue
                else:#have not seen the data yet
                    reactiontree.strand_hash_table[datastrand]={}
                    reactiontree.strand_hash_table[datastrand][strandIDmod]=1
                    reactiontree.inc_node_count()
    #create a unoptimized tree lite to calculate necessary padding overheads
    unopt_tree=construct_tree_unoptimized_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy,h_array)
    reactiontree.add_pad_count(unopt_tree.get_pad_count()/len(strands)) #can optimize therefore the price of padding count should be only paid once, padding is also not data dependent only alignment dependent per strand
    reactiontree.dec_node_count(unopt_tree.get_singleton_nodes_removed()) #account for removal of nodes
    return reactiontree


def construct_tree_transform_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None):
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    for strand_index, s in enumerate(strands):
        #print len(s)
        h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
        skip_array=[False]*strand_length_in_codewords #If an entry is True: you skip this sub strand
        for height in range(1,h+1)[::-1]:
            block_group_size=(reactiontree.num_overhangs-1)**height
            for strandID, i in enumerate(range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length))):#breaks up the input strand into sub reactions
                strandIDmod=strandID%num_overhangs
                if i == (len(s)-(reactiontree.overhang_length)): break
                if skip_array[strandID*block_group_size]: continue #this tree is already optimized out
                end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
                if end_index>=len(s): end_index=len(s)
                substrand=s[i:end_index]
                datastrand=""
                #cut data out of substrand
                for data_start in range(reactiontree.overhang_length,len(substrand),reactiontree.codeword_length+reactiontree.overhang_length):
                    datastrand+=substrand[data_start:data_start+reactiontree.codeword_length]
                #print len(datastrand)
                #print codeword_length*block_group_size
                assert len(datastrand)==codeword_length*block_group_size
                if datastrand in reactiontree.strand_hash_table:
                    if strandIDmod in reactiontree.strand_hash_table[datastrand]:#check overhang requirement
                        #have a match, take the height statistic, don't count the reaction, mark the skip_array
                        if reactiontree.h_array is not None:
                            reactiontree.h_array[height-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
                        for _ in range(0,block_group_size):
                            skip_array[_+strandID*block_group_size]=True #set the underlying codewords to all True
                        if height>1 and num_overhangs>=5:
                            #check to see if this is a product of a transform
                            transform_info=reactiontree.strand_hash_table[datastrand][strandIDmod]
                            if transform_info[0]>0:#are using the product of a transform
                                if transform_info[1]==0:
                                    reactiontree.inc_node_count() #activate transform
                                    if transform_info[0]==2:
                                        transform_info_2=reactiontree.strand_hash_table[datastrand][transform_info[2]]
                                        if transform_info_2[0]>0:
                                            if transform_info_2[1]==0:
                                                reactiontree.inc_node_count()
                                                reactiontree.strand_hash_table[datastrand][transform_info[2]][1]+=1
                            reactiontree.strand_hash_table[datastrand][strandIDmod][1]+=1 #increment use count
                    
                        elif height==1 or num_overhangs<5:
                            reactiontree.strand_hash_table[datastrand][strandIDmod][1]+=1 #increment use count     
                        continue #move onto next substrand

                    else: #have seen the data, but not this strandIDmod (overhang version)
                        if num_overhangs<5 or height==1:
                            reactiontree.strand_hash_table[datastrand][strandIDmod]=[0,1]
                            reactiontree.inc_node_count()
                        elif num_overhangs>=5 and height>1: #should have transforms upon first creation of the data, assert
                            assert 0
                        continue
                else:#have not seen the data yet
                    reactiontree.strand_hash_table[datastrand]={}
                    reactiontree.strand_hash_table[datastrand][strandIDmod]=[0,1,0]
                    reactiontree.inc_node_count()
                    if height>1 and num_overhangs>=5: #tentatively create the necessary transforms for all possible strandIDmods using the current strandIDmod as the seed
                        if strandIDmod==0:
                            strandIDmod_b=num_overhangs-1
                        else:
                            strandIDmod_b=strandIDmod-1
                        strandIDmod_a=(strandIDmod+1)%num_overhangs
                        strandIDstart=strandIDmod_a
                        while strandIDstart!=strandIDmod:
                            if strandIDstart==strandIDmod_a or strandIDstart==strandIDmod_b:
                                transform_info=[2,0,(strandIDmod+3)%num_overhangs] #tuple of [depth of transform, and trasform use count]
                                reactiontree.strand_hash_table[datastrand][strandIDstart]=transform_info
                            else:
                                transform_info=[1,0]
                                reactiontree.strand_hash_table[datastrand][strandIDstart]=transform_info
                            strandIDstart=(strandIDstart+1)%num_overhangs

    #create a unoptimized tree lite to calculate necessary padding overheads
    unopt_tree=construct_tree_unoptimized_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy,h_array)
    reactiontree.add_pad_count(unopt_tree.get_pad_count()/len(strands)) #can optimize therefore the price of padding count should be only paid once, padding is also not data dependent only alignment dependent per strand
    reactiontree.dec_node_count(unopt_tree.get_singleton_nodes_removed()) #account for removal of nodes
    return reactiontree



#analysis for inserting rotation reactions in order to ultimately reduce reactions, opt_dictionary is that of 
def construct_tree_rotate_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy=None,h_array=None,opt_dictionary=None):
    DIV=4
    reactiontree=ReactionTree(num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,repair_strategy)
    total_NOP_reactions=0
    total_strand_NOP_inserts=0 #tracks the number of NOP inserts in the strand
    for strand_index, s in enumerate(strands):
        h=int(math.ceil(math.log(strand_length_in_codewords,num_overhangs-1)))
        skip_array=[False]*strand_length_in_codewords #If an entry is True: you skip this sub strand
        shift_array=[0]*strand_length_in_codewords #indicates the number that should be added to strandIDmod in order to derive new strandID
        insert_height_array=[0]*(h+200) #this array contains the number of NOP reactions inserted at a certain height
        strand_NOP_inserts=0 #tracks the number of NOP inserts in the strand
        for height in range(1,h+1)[::-1]:
            block_group_size=(reactiontree.num_overhangs-1)**height
            NOP_inserts=0
            for strandID, i in enumerate(range(0,len(s),block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length))):#breaks up the input strand into sub reactions
                if i == (len(s)-(reactiontree.overhang_length)): break
                if skip_array[strandID*block_group_size]: continue #this tree is already optimized out
                opt_strandID=strandID%num_overhangs #this strand ID preserves its original ID 
                strandIDmod=(shift_array[strandID*block_group_size]+strandID)%num_overhangs
                end_index=i+((block_group_size*(reactiontree.codeword_length+reactiontree.overhang_length)+reactiontree.overhang_length))
                if end_index>=len(s): end_index=len(s)
                substrand=s[i:end_index]
                datastrand=""
                #cut data out of substrand
                for data_start in range(reactiontree.overhang_length,len(substrand),reactiontree.codeword_length+reactiontree.overhang_length):
                    datastrand+=substrand[data_start:data_start+reactiontree.codeword_length]
                assert len(datastrand)==codeword_length*block_group_size
                opt_use_count=opt_dictionary[datastrand][opt_strandID]
                if datastrand in reactiontree.strand_hash_table:
                    if strandIDmod in reactiontree.strand_hash_table[datastrand]:#check overhang requirement
                        #have a match, take the height statistic, don't count the reaction, mark the skip_array
                        if reactiontree.h_array is not None:
                            reactiontree.h_array[height-1]+=1 #increment counter tracking the number of times redundancy was found at a certain height
                        for _ in range(0,block_group_size):
                            skip_array[_+strandID*block_group_size]=True #set the underlying codewords to all True
                        continue #move onto next substrand
                    if opt_use_count>10:
                        #print "too large of use"
                        #this has been seen to be used a lot just make the node
                        reactiontree.strand_hash_table[datastrand][strandIDmod]=True
                        reactiontree.inc_node_count()                        
                    else: #have seen the data, but not this strandIDmod (overhang version)
                        #try to rotate to the current strandIDmod
                        #print "going to rotate !!!!!!!!!!!!!!!!!!"
                        strandIDmod_e=[_ for _ in reactiontree.strand_hash_table[datastrand]]
                        distance=float('inf')
                        for dest in strandIDmod_e: #find smallest distance
                            #get the distance (source to destination distance)  in which we are rotating
                            if strandIDmod < dest:
                                _distance=dest-strandIDmod
                            else:
                                assert strandIDmod!=dest
                                distance_to_0=num_overhangs-strandIDmod
                                distance_from_0=dest
                                _distance=distance_to_0+distance_from_0
                            assert _distance>0
                            if _distance<distance:
                                distance=_distance
                        distance=int(distance)
                        if distance < (num_overhangs-1)/DIV:  #minimize distance to lower overhead impact
                            #print " strand index {} height {} distance {} num overhangs: {} num_ovrehangs-1/3 {}".format(strand_index,height,distance,num_overhangs,(num_overhangs-1)/3)
                            #increment the insertion counters
                            for _ in range(0,len(shift_array[block_group_size*strandID::])):
                                shift_array[_+strandID*block_group_size]+=distance                                
                            NOP_inserts+=distance #track NOP inserts at this height
                        else:
                            reactiontree.strand_hash_table[datastrand][strandIDmod]=True
                            reactiontree.inc_node_count()
                    continue
                else:#have not seen the data yet
                    reactiontree.strand_hash_table[datastrand]={}
                    reactiontree.strand_hash_table[datastrand][strandIDmod]=True
                    reactiontree.inc_node_count()
        
            insert_height_array[height-1]+=NOP_inserts
            for _index, _ in enumerate(shift_array):
                shift_array[_index]=shift_array[_index]*(num_overhangs-1)#increase the shift amount to be consistant down the tree
            strand_NOP_inserts+=NOP_inserts
        total_strand_NOP_inserts+=strand_NOP_inserts
        #add analysis of reactions to add for this strand
        height_array_index=0
        #print insert_height_array
        #print height_array_index
        while insert_height_array[height_array_index]==0:
            height_array_index+=1
            if height_array_index==len(insert_height_array): break
        if height_array_index==len(insert_height_array): continue 
        #have the start of added rotation reactions
        reaction_counter=0
        #print insert_height_array
        while True:
            #keep tracking new reactions up the tree until the added reaction at the next level is
            local_reaction_counter=math.ceil(float(insert_height_array[height_array_index])/float((num_overhangs-1)))
            reaction_counter+=local_reaction_counter
            insert_height_array[height_array_index+1]+=local_reaction_counter
            if insert_height_array[height_array_index+1]-local_reaction_counter==0 and local_reaction_counter==1:#adding 1 to next layer that had nothing 
                break
            height_array_index+=1#counter for adding in last reaction to merge all NOPS
        reaction_counter+=1
        total_NOP_reactions+=reaction_counter
    print("Total NOP Added {} Total Strands {} Avg NOPs/strand {}".format(total_strand_NOP_inserts,len(strands),float(total_strand_NOP_inserts)/float(len(strands))))
    reactiontree.add_node_count(total_NOP_reactions) #add back in NOP counts insertedd      
    #create a unoptimized tree lite to calculate necessary padding overheads
    unopt_tree=construct_tree_unoptimized_lite(strands,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,repair_strategy,h_array)
    reactiontree.add_pad_count(unopt_tree.get_pad_count()/len(strands)) #can optimize therefore the price of padding count should be only paid once, padding is also not data dependent only alignment dependent per strand
    reactiontree.dec_node_count(unopt_tree.get_singleton_nodes_removed()) #account for removal of nodes
    return reactiontree











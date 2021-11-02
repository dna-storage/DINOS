'''
Author: Kevin Volkel

Filename: base_node.py

Description: This file contains the base class that represents a node in an overhang assembly graph

'''
import math
import overhang.dnastorage_utils.codec.base_conversion as bc #support for base conversion, needed to initialize lookup table


class ReactionNode:#this class acts as a container for information relevant to nodes in a overhang assembly graph
    def __init__(self, strandID, height,is_pad):

        self._strandIDList=[] #captures the strand IDs that use this node
        if type(strandID) is not list:
            self._strandIDList.append(strandID)
        else:
            self._strandIDList+=strandID #inputting a list
        self._use_count=1 #indicates the number of users of this node
        self._sub_nodes=[] #list of children subnodes, if node is a base node, set this to None
        self._visited=0
        self._is_pad=is_pad #inidcates wheter or not the node is used for padding
        self._is_terminator=0
        self._number_true_children=0 #number of true children (non padding) a node has
        self._height=height #height in the tree that the node is, nodes should only be shared amongst other nodes at the same height, otherwise implies sharing an incomplete reaction (may be ok to share the null padded strand)
        self._pad_strand="" #padding added to the node to make sure overhangs are correct
        self._corrected_substrand=""
    #methods to update node characteristics upon creation
    def inc_use_count(self):
        self._use_count+=1
    def dec_use_count(self):
        self._use_count-=1
    def get_use_count(self):
        return self._use_count

    def get_corrected_substrand(self):
        return self._corrected_substrand
    def set_corrected_substrand(self,corrected_strand):
        self._corrected_substrand=corrected_strand
    
    #methods to manage what strand indexes are using this reaction
    def insert_strand_ID(self,ID):
        if type(ID) is not list:
            self._strandIDList.append(ID)
        else:
            self._strandIDList+=ID

    def is_ID(self,ID):
        if ID in self._strandIDList:
            return True
        else:
            return False

    #methods for toggling the visited state of the node
    def set_visited(self):
        self._visited=1
    def clear_visited(self):
        self._visited=0
    def get_visited(self):
        return self._visited
    def invert_visited(self):
        if self._visited==0:
            self._visited=1
        elif self._visited==1:
            self._visited=0

    #methods to manage node children
    def add_child(self,child):#add child node to this node
        self._sub_nodes.append(child)
    def get_child_list(self):
        return self._sub_nodes

    #methods to manage true children
    def set_true_children(self,number_true):
        self._number_true_children=number_true
    def get_true_children(self):
        return self._number_true_children
    #methods to manage terminator status
    def set_term(self):
        self._is_terminator=1
    def is_term(self):
        return self._is_terminator
    def clear_term(self):
        self._is_terminator=0

    #methods to manage whether a node is padding or not
    def set_pad(self):
        self._is_pad=1
    def clear_pad(self):
        self._is_pad=0
    def is_pad(self):
        return self._is_pad

    #methods to manage padding strand 
    def set_pad_strand(self,pad_strand):
        self._pad_strand=pad_strand
    def get_pad_strand(self):
        return self._pad_strand
    

    #methods to manage height information
    def set_height(self,height):
        self._height=height
    def get_height(self):
        return self._height


class analytical_ReactionTree: #analytics based reaction tree so that we don't use so much memory in non-optimization cases
    def __init__(self,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,end_repair):
        self._tree_stats={}
        #internal state for tree walking status and stats
        self._tree_stats["insert_per_strand"]=0
        self._tree_stats["node_count"]=0
        self._tree_stats["pad_node_count"]=0 #nodes needed to pad tree for overhang correctness
        self._tree_stats["singleton_nodes_removed"]=0 #number of singeton nodes removed
        
        #state describing the properties of the reaction network being constructed
        #self.strand_length=strand_length
        self.strand_length_in_codewords=strand_length_in_codewords
        self.codeword_length=codeword_length
        self.overhang_length=overhang_length
        self.ideal_children_nodes=num_overhangs-1 #number of nodes that a child should have
        self.end_repair_technique=end_repair # end repair technique used, 0: direct end repair (m^2) extra strands, 1: same reaction repair (2m extra strands) 2: reaction repair (m extra strands)
        self.strand_hash_table={} #hash table used to track redundancy
        self.num_overhangs=num_overhangs
        self.overhang_length=int(math.ceil(math.log(num_overhangs,4))) #overhangs are represented as base 4 numbers
    
    def add_node_count(self,add):
        self._tree_stats["node_count"]+=add #can increment manually to avoid walking
        
    def minus_node_count(self, dec_count):#decrement the count of nodes 
        self._tree_stats["node_count"]-=dec_count

    def add_pad_count(self,add):
        self._tree_stats["pad_node_count"]+=add

    def minus_pad_count(self,dec_count):#decrement the count of padding nodes
        self._tree_stats["pad_node_count"]-=dec_count

    def set_insert_per_strand(self,insert_count):
        self._tree_stats["insert_per_strand"]=insert_count
        
    def order(self):
        return self._tree_stats["node_count"]

    def get_pad_count(self):
        return self._tree_stats["pad_node_count"]

    def set_singleton_nodes_removed(self,num):
        self._tree_stats["singleton_nodes_removed"]=num

    def get_singleton_nodes_removed(self):
        return self._tree_stats["singleton_nodes_removed"]
    
    
class ReactionTree: #this class will take the place of using the networkx module, should be more lightweight and faster
    
    def __init__(self,num_overhangs,overhang_length,codeword_length,strand_length_in_codewords,h_array,end_repair):
        self.num_to_overhang=[] #array used to lookup the string representing an overhang using a number, should be length m.
        self.overhang_to_num={} #find the ID from the overhang

        #internal state for tree walking status and stats
        self._term_nodes=[] #list of terminating nodes
        self._tree_stats={}
        self._tree_stats["node_count"]=0
        self._tree_stats["pad_node_count"]=0 #nodes needed to pad tree for overhang correctness
        self._tree_walk_status=0
        self._tree_walk_count=0

        #state describing the properties of the reaction network being constructed
        #self.strand_length=strand_length
        self.strand_length_in_codewords=strand_length_in_codewords
        self.codeword_length=codeword_length
        self.overhang_length=overhang_length
        self.ideal_children_nodes=num_overhangs-1 #number of nodes that a child should have
        self.end_repair_technique=end_repair # end repair technique used, 0: direct end repair (m^2) extra strands, 1: same reaction repair (2m extra strands) 2: reaction repair (m extra strands)
        self.strand_hash_table={} #hash table used to track redundancy
        self.h_array=h_array #stat array used to track at what height data is shared
        self.num_overhangs=num_overhangs

        self.num_to_overhang=[]
        self.overhang_length=int(math.ceil(math.log(num_overhangs,4))) #overhangs are represented as base 4 numbers
        
        #initilaize overhang lookup tables
        self._initialize_lookup_table(num_overhangs)
        
    def _initialize_lookup_table(self,num_overhangs):
        self.num_to_overhang=[bc.convertQuarnary(_,self.overhang_length)[::-1] for _ in range(0,self.num_overhangs)] #arbitrary strings that ID overhangs
        for overhangID,overhang in enumerate(self.num_to_overhang):
            self.overhang_to_num[overhang]=overhangID

    
    def tree_walk(self,node=None): #walks the tree recursively, figuring out node stats
        if node==None:#top of the tree
            for strand_terminator in self._term_nodes:
                if (strand_terminator.get_visited()==0 and self._tree_walk_status==0) or (strand_terminator.get_visited()==1 and self._tree_walk_status==1) :
                    tree_walk(node=strand_terminator)
                    strand_terminator.invert_visited() #do inverse of node visited, that way we don't have to clear the visited check on a subsequent walk, we can just walk the tree
                    self._tree_stats["node_count"]+=1
        else:
            for child in node.get_child_list():
                if (child.get_visited()==0 and self._tree_walk_status==0) or (strand_terminator.get_visited()==1 and self._tree_walk_status==1):
                    tree_walk(node=child)
                    child.invert_visited()
                    self._tree_stats["node_count"]+=1
                    
    def clear_stats(self):#clear stats
        for stat in self._tree_stats:
            self._tree_stats[stat]=0
            
    def collect_stats(self):
        self.clear_stats()
        self.tree_walk()
        #make sure the status is inverted
        if self._tree_walk_status==0:
            self._tree_walk_status=1
        else:
            self._tree_walk_status=0
        self._tree_walk_count+=1

    def order(self):
        assert (self._tree_walk_count>0 or self._tree_stats["node_count"]>0) #make sure we walk the tree at least once
        return self._tree_stats["node_count"]
    

    def add_terminator(self,node):
        self._term_nodes.append(node)

    def inc_node_count(self):
        self._tree_stats["node_count"]+=1 #can increment manually to avoid walking

    def add_node_count(self,add):
        self._tree_stats["node_count"]+=add
    
    def dec_node_count(self, dec_count):#decrement the count of nodes 
        self._tree_stats["node_count"]-=dec_count

    def inc_pad_count(self):
        self._tree_stats["pad_node_count"]+=1

    def add_pad_count(self,add):
        self._tree_stats["pad_node_count"]+=add

    def dec_pad_count(self,dec_count):#decrement the count of padding nodes
        self._tree_stats["pad_node_count"]-=dec_count

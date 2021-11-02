from base_conversion import *
from base import *
import math

class BinaryStringCodec(BaseCodec): #simple class used to create binary string code word arrays 
    def __init__(self,numberBytes,bits_per_block,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj,Policy=Policy)
        self._numberBytes=numberBytes
        self._bits_per_block=bits_per_block
    def _encode(self,strand):
        enc_strand=[ convertBytetoBinary(c,int(self._bits_per_block),reverse=True) for c in strand] #should get an array of binary strings, each with length bits_per_block
        #print enc_strand
        '''
        #convert enc_strand to a list of binary cahracters
        binary_string=[]
        for b in enc_strand:
            binary_string+=list(b)
        codeword_list=[]
        #consolidate bits into codewords indicated by bits_per_block
        for binary_strart_index in range(0,len(binary_string),self._bits_per_block):
            codeword_list.append("".join(binary_string[binary_strart_index:binary_strart_index+self._bits_per_block]))
        return codeword_list
        '''
        return enc_strand #already have the codewords at the right length, no longer need to cut up binary string

class InsertOverhangs(BaseCodec): #this class creates a set of arbitrary overhangs to be iserted between codewords
    def __init__(self,num_overhangs,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj,Policy=Policy)
        self._num_overhangs=num_overhangs
        self.overhang_length=int(math.ceil(math.log(self._num_overhangs,4)))
        self._overhangs=[convertQuarnary(_,self.overhang_length)[::-1] for _ in range(0,self._num_overhangs)] #arbitrary strings that ID overhangs
    def _encode(self,strand): #insert between each block
        overhang_counter=1
        overhang_strand=[]
        overhang_strand.append(self._overhangs[0])
        for index, building_block in enumerate(strand):
            overhang_strand.append(building_block)
            overhang_strand.append(self._overhangs[overhang_counter])
            overhang_counter+=1
            overhang_counter%=self._num_overhangs #rotate overhang counter
        return overhang_strand #should have building blocks with overhangs inserted in between

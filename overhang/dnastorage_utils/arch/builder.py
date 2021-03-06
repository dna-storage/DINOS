from dnastorage.codec.phys import *
from dnastorage.codec.strand import *
from dnastorage.codec.block import *
from dnastorage.codec.LayeredCodec import *
from dnastorage.exceptions import *
from dnastorage.codec.binarystringcodec import *
import math 
def available_file_architectures():
    return ['OverhangBitStringStrand']



def build_overhang_bitstring_strand(is_enc,pf,primer5,primer3,strand_length,num_overhangs,bits_per_block): #creates strands of length strand_length which is measured in codewords, bits_per_block controls the size of a codeword, and the number of overhangs indicates the number of connecting regions available
    
    pol=AllowAll()
    if is_enc==True:
        #block indexes meausred in number of blocks
        interBlockIndex=int(math.ceil(float(8*3)/float(bits_per_block))) #24 bits worth of strand indexing
        intraBlockIndex=int(math.ceil(float(8*0)/float(bits_per_block))) #0 bits of block indexing
        index=intraBlockIndex+interBlockIndex
        #This is a simple implementation with no error correction built in, simply reads in a strandlength worth of data and then chooses the sequence of overhangs to use
        #resultant strand from the encoder will follow the format Overhang:Data:Overhang ...

        #blank outer encoding, no ECC placed on outer encoding 
        blockCodec=DoNothingOuterCodec(packetSize=strand_length,
                                       payloadSize=strand_length,
                                       Policy=pol)


        blockToStrand = BlockToStrand(strand_length,strand_length,Policy=pol,\
                                  intraIndexSize=intraBlockIndex,\
                                  interIndexSize=interBlockIndex)

        blockToStrand.codewordSize=bits_per_block #tell the block to strand layer the size of codewords we are using, important for index conversions
        
        strandCodec=ReedSolomonInnerCodec(0,Policy=pol)#inner ecc set to 0, we need to make a long ECC for inner strand protection of long strands

        codewords=BinaryStringCodec(strand_length+index,bits_per_block,Policy=pol)#converts to a binary string

        overhangs=InsertOverhangs(num_overhangs,Policy=pol,CodecObj=codewords) #overhangs inserts overhang sequences in between codewords
        
        pre = PrependSequence(primer5,isPrimer=True,Policy=pol)
        
        app = AppendSequence(reverse_complement(primer3), CodecObj=pre, isPrimer=True)
        
        physCodec= app #append and prepend primer strings
        
        pf.codewordSize=bits_per_block #set the packetized file codeword size
        
        enc=LayeredEncoder(pf,blockSizeInBytes=strand_length,strandSizeInBytes=strand_length,\
                           blockCodec=blockCodec,\
                           blockToStrandCodec=blockToStrand,\
                           strandCodec=strandCodec,\
                           strandToCodewordCodec=overhangs,\
                           codewordToPhysCodec=CombineCodewords(),\
                           physCodec=physCodec,Policy=pol)
        return enc #result will be a a set of binary strings with overhang tags in between


def build_encode_architecture(arch, pf, primer5, primer3,num_overhangs=None,codeword_size=8):
    if arch == 'OverhangBitStringStrand':
        return build_overhang_bitstring_strand(True,pf,primer5,primer3,strand_length=1000,num_overhangs=num_overhangs,bits_per_block=codeword_size)#encoding algorithm creates an encoding file with will be interpreted by a compiler to create instructions that will create the strand 
    

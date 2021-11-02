from dnastorage.codec.phys import CombineCodewords
from dnastorage.codec.block import *
from dnastorage.codec.codecfile import *
from math import ceil,log

class LayeredEncoder(EncodePacketizedFile):

    def __init__(self,packetizedFile,minIndex=0,\
                 strandSizeInBytes=10,blockSizeInBytes=200*10,\
                 blockCodec=None,\
                 blockToStrandCodec=None,\
                 strandCodec=None,\
                 strandToCodewordCodec=None,\
                 codewordToPhysCodec=CombineCodewords(),\
                 physCodec=None,\
                 Policy=None):
        
        EncodePacketizedFile.__init__(self,packetizedFile,minIndex=minIndex)
        # set packetSize (can't hurt to do it again here)
        self._packetizedFile.packetSize = blockSizeInBytes

        # make sure sizes are consistent with each other
        # blocks must be whole multiple of strand size
        assert blockSizeInBytes % strandSizeInBytes == 0

        self.blockCodec = blockCodec
        self.blockToStrandCodec = blockToStrandCodec
        self.strandCodec = strandCodec
        self.strandToCodewordCodec = strandToCodewordCodec
        self.codewordToPhysCodec = codewordToPhysCodec
        self.physCodec = physCodec
        self.strandSizeInBytes = strandSizeInBytes
        self.blockSizeInBytes = blockSizeInBytes
        return

    def _layered_encode(self, block):
        # perform block encoding
        enc_block = self.blockCodec.encode(block)

        # convert block into strands
        strands = self.blockToStrandCodec.encode(enc_block)


        # protected logical strands
        final_strands = []
        tmp_strands = []
        for s in strands:
            ecc_s = self.strandCodec.encode(s)
            tmp_strands.append(ecc_s)
            cw_s = self.strandToCodewordCodec.encode(ecc_s)
            # phys codecs expect a DNA sequnce as a string
            phys_s = self.codewordToPhysCodec.encode(cw_s)            
            final = self.physCodec.encode(phys_s)
            final_strands.append(final)

        return final_strands

    
    def dummy_encode(self):
        """ Dummy encode let's us run the sequence of codecs and produce a single
            strand. We can look at that strand to determine its length or possibly
            other properties. Helpful when designing codecs.
        """
        block = [ randint(0,255) for _ in range(self.blockSizeInBytes) ]
        strands = self._layered_encode( (randint(0,255**2),block) )
        return strands
                
    def encode(self):        
        block = self._encode()
        try:
            block = (block[0],[ ord(_) for _ in block[1] ])
        except Exception as e:
            # This is pretty horrible, kids don't try this at home!
            block = (block[0],[ _ for _ in block[1] ])
        return self._layered_encode(block) # get entire block
        

class LayeredDecoder(DecodePacketizedFile):
    
    def __init__(self,packetizedFile,minIndex=0,\
                 strandSizeInBytes=10,blockSizeInBytes=200*10,\
                 blockIndexSize=2,
                 physCodec=None,\
                 physToStrandCodec=None,\
                 strandCodec=None,\
                 strandToBlockCodec=None,\
                 blockCodec=None,\
                 Policy=None,
                 intraBlockIndexSize=1):\
                         
        DecodePacketizedFile.__init__(self,packetizedFile,minIndex=minIndex)
        # set packetSize (can't hurt to do it again here)
        self._packetizedFile.packetSize = blockSizeInBytes

        # make sure sizes are consistent with each other
        # blocks must be whole multiple of strand size
        assert blockSizeInBytes % strandSizeInBytes == 0

        self.physCodec = physCodec
        self.physToStrandCodec = physToStrandCodec
        self.strandCodec = strandCodec
        self.strandToBlockCodec = strandToBlockCodec        
        self.blockCodec = blockCodec
                
        self.strandSizeInBytes = strandSizeInBytes
        self.blockSizeInBytes = blockSizeInBytes
        self.all_strands = []
        self.blockIndexSize = blockIndexSize
        self.intraBlockIndexSize = intraBlockIndexSize
        self._Policy = Policy
        return


    def _layered_decode_phys_to_strand(self, phys_strand):
        # perform block encoding
        try:
            phys_s = self.physCodec.decode(phys_strand)
            cw_s = self.physToStrandCodec.decode(phys_s)
            s = self.strandCodec.decode(cw_s)
        except DNAStorageError as p:
            s = [-1] + [ 0 for _ in range(self.strandSizeInBytes-1) ]
        
        return s

    def decode_from_phys_to_strand(self, s):
        return self._layered_decode_phys_to_strand(s)
    
    def decode(self, phys_strand, bypass=False, input_key=None, input_value=None):
        if bypass==False:
            try:            
                s = self._layered_decode_phys_to_strand(phys_strand)
                self.all_strands.append(s)
            except DNAMissingPrimer as p:
                # just ignore strands that don't have a required primer
                pass
            except DNAStorageError as e:
                if self._Policy.allow(e):
                    pass
                else:
                    raise e
        else:
            self.all_strands.append(input_value)

    def _attempt_final_decoding(self):
        # do voting here!!
        self.voted_strands = doMajorityVote(self.all_strands,\
                                            self.blockIndexSize+self.intraBlockIndexSize)
        if self.blockIndexSize>0:
            blocks = partitionStrandsIntoBlocks(self.all_strands,self.blockIndexSize)
        else:
            blocks = [ (0,self.all_strands) ]
        #print "number of blocks=",len(blocks)
        for b in blocks:
            idx = b[0]
            try:
                #print "attempt",idx,len(b[1])
                b_contig = self.strandToBlockCodec.decode(b)
                #print "attempt",b_contig[0],len(b_contig[1])
                b_noecc = self.blockCodec.decode(b_contig)
                #print "attempt",b_noecc[0],len(b_noecc[1])
                self.writeToFile(b_noecc[0],b_noecc[1])

            except DNAStorageError as e:
                if self._Policy.allow(e):
                    continue                
                else:
                    raise e
            except Exception as e:
                print (str(e))
                print (b)
                
            
    def write(self):
        self._attempt_final_decoding()
        super(LayeredDecoder,self).write()
            
    def dummy_write(self):
        self._attempt_final_decoding()
        return self._packetizedFile.dummy_write()

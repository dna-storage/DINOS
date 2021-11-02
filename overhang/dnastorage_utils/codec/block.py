from dnastorage.codec.base import *
from math import log, ceil
from dnastorage.codec import base_conversion
from dnastorage.codec.reedsolomon.rs import ReedSolomon,get_reed_solomon,ReedSolomonError
from collections import Counter

import logging
logger = logging.getLogger('dna.storage.codec.block')
logger.addHandler(logging.NullHandler())


class BlockToStrand(BaseCodec):
    """ This class encodes an (index,block) pair into a strand that contains 
        it's own index follwed by the data.  After this conversion, the index
        simpy appears as part of the strand.
        The _encode returns a list of strands.
    
        The _decode function expects to receive a set of strands that belong to
        the same block. Then they are assembled into a cohesive block and returned.     
    """

    def __init__(self, strandSizeInBytes, blockSizeInBytes, intraIndexSize=1, interIndexSize=2, nSyms=256, CodecObj=None, Policy=None, filterZeroes=False):
        super(BlockToStrand,self).__init__(CodecObj=CodecObj,Policy=Policy)
        self._strandSizeInBytes = strandSizeInBytes 
        self._blockSize = blockSizeInBytes
        self._interIndex = interIndexSize
        self._intraIndex = intraIndexSize
        self._nSyms = nSyms
        self.__codewordSize=8 #default codeword size to a byte to support legacy code
        self._removeZeroes = filterZeroes
        logger.info("filterZeroes = {}".format(filterZeroes))



    def __set_codewordSize(self,val):
        self.__codewordSize=val
    def __get_codewordSize(self):
        return self.__codewordSize
    codewordSize=property(__get_codewordSize,__set_codewordSize)
    
    def allzeroes(self, l):
        for z in l:
            if z != 0:
                return False
        return True

    def _filter_zeroes(self, block):
        if self._removeZeroes == False:
            return [ False for _ in range(0,len(block),self._strandSizeInBytes) ]

        zeroes = []
        for i in range(0,len(block),self._strandSizeInBytes): 
            strand = block[i:i+self._strandSizeInBytes]
            if self.allzeroes(strand):
                zeroes.append(True)
            else:
                zeroes.append(False)
                    
        fil = [ False for _ in range(len(zeroes)) ]
        assert len(fil) == len(zeroes)
        for i,z in enumerate(zeroes):
            if i > 0 and i < len(zeroes)-1:
                if zeroes[i-1]==True and zeroes[i+1]==True and zeroes[i]==True:
                    fil[i] = True

        #print fil
        return fil
                    
    def _encode(self, packet):
        strands = []
        #bindex = base_conversion.convertIntToBytes(packet[0],self._interIndex)
        #make the base conversion of int into arbitrary codeword size
        bindex=base_conversion.convertIntToCodeword(packet[0],self._interIndex,codewordSize=self.codewordSize)
        block = packet[1]
        #print block

        fil = self._filter_zeroes(block)
            
        #assert len(block) <= self._blockSize
        
        assert len(bindex) <= self._interIndex
        
        for i in range(0,len(block),self._strandSizeInBytes):
            #if allzeroes(block[i:i+self._strandSizeInBytes]):
            #    continue
            if fil[int(i/self._strandSizeInBytes)]==True: # don't emit zeroes
                stats.inc("BlockToStrand::filterAllZeroStrand")
                continue
            
            sindex = base_conversion.convertIntToCodeword(i/self._strandSizeInBytes,self._intraIndex,codewordSize=self.codewordSize)
            assert len(sindex) <= self._intraIndex
            #print sindex
            s = bindex + sindex  + block[i:i+self._strandSizeInBytes]            
            strands.append(s)
        return strands

    def _decode(self, packet):        
        """ accepts a list of strands and coverts them into a contiguous block """
        d = {} # track all strands we find
        
        prior_bindex = packet[0]
        strands = packet[1]
        max_index = 0
        for s in strands:
            bindex = base_conversion.convertBytesToInt(s[0:self._interIndex])
            sindex = base_conversion.convertBytesToInt(s[self._interIndex:self._interIndex+self._intraIndex])

            if bindex != prior_bindex:
                e = DNABlockBadIndex("Bad index {}!={}".format(bindex,prior_bindex))
                #print "Bad index {}!={}".format(bindex,prior_bindex)
                if self._Policy.allow(e):
                    # ignore this strand
                    if d.has_key(sindex):
                        # already found a strand with this sindex, don't use this one
                        continue
                    else:
                        pass
                else:
                    raise e

            if len(s[self._interIndex+self._intraIndex:]) != self._strandSizeInBytes:            
                # handle exception condition here!
                err = DNABlockPayloadWrongSize("Strand payload size should be {} bytes but was {}.".format(self._strandSizeInBytes,len(s[self._interIndex+self._intraIndex:])))
                if self._Policy.allow(err):
                    if len(s[self._interIndex+self._intraIndex:]) > self._strandSizeInBytes:
                        d[sindex] = s[self._interIndex+self._intraIndex:self._interIndex+self._intraIndex+self._strandSizeInBytes]
                    else:
                        d[sindex] = s[self._interIndex+self._intraIndex:]+[-1]*(self._strandSizeInBytes-len(s[self._interIndex+self._intraIndex:]))
                else:
                    raise err
            else:
                d[sindex] = s[self._interIndex+self._intraIndex:]
                
        data = []
        max_index = self._blockSize / self._strandSizeInBytes
        for i in range(max_index):
            if not d.has_key(i):
                if self._removeZeroes:
                    if not (i-1 in d) and not ((i+1) in d):
                        # guess that it's all zeros
                        stats.inc("BlockToStrand::decode::guessAllZeroStrand")
                        data += [ 0 for _ in range(self._strandSizeInBytes) ]
                    elif (i-1) in d and self.allzeroes(d[i-1]) and not (i+1 in d):
                        # guess that it's all zeros
                        stats.inc("BlockToStrand::decode::guessAllZeroStrand")
                        data += [ 0 for _ in range(self._strandSizeInBytes) ]
                    elif (i+1) in d and self.allzeroes(d[i+1]) and not (i-1 in d):
                        # guess that it's all zeros
                        stats.inc("BlockToStrand::decode::guessAllZeroStrand")
                        data += [ 0 for _ in range(self._strandSizeInBytes) ]
                    else:
                        err = DNABlockMissingIndex("Block missing index = {}".format(i))
                        stats.inc("BlockToStrand::decode::guessMissingStrand")
                        if self._Policy.allow(err):
                            #print "Block missing index = {}".format(i)
                            data += d.get(i,[-1 for _ in range(self._strandSizeInBytes)])
                        else:
                            raise err
                else:
                    err = DNABlockMissingIndex("Block missing index = {}".format(i))
                    if self._Policy.allow(err):
                        #print "Block missing index = {}".format(i)
                        data += d.get(i,[-1 for _ in range(self._strandSizeInBytes)])
                    else:
                        raise err
            else:
                data += d[i]
                
        # smaller is allowed due to end of file, larger is a potential problem
        if len(data) > self._blockSize:
            err = DNABlockTooLargeError("Block too large, is {} should be {}".format(len(data),self._blockSize))
            if self._Policy.allow(err):                
                data = data[:self._blockSize]
            else:
                raise err

        #print data
        return prior_bindex,data
        

class DoNothingOuterCodec(BaseCodec):
    """This is a simple do nothing Outer Codec to bypass this step when testing new designs"""
    def __init__(self,packetSize,payloadSize,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._packetSize=packetSize
        self._payloadSize=payloadSize
    def _encode(self,packet):
        """simply check length of packet to make sure it is a multiple of payload size"""
        index = packet[0]
        data = packet[1][:]
        
        assert len(packet[1])<=self._packetSize

        # hopefully the last packet if it's not a full packetSize
        # so, pad it out with zeros to make an even number of strands if
        # isn't already
        
        if len(data) < self._packetSize:
            stats.inc("DoNothingOuterCodec.padPacket")
            # normalize to multiple of payloadSize
            rem = len(data) % self._payloadSize
            data += [0]*(self._payloadSize-rem)

        #print "Data length {}".format(len(data))
        assert len(data)==self._packetSize
        return packet[0],data

        
class ReedSolomonOuterCodec(BaseCodec):
    """ReedSolomonOuterCodec takes a block of data as input and produces a new
    block of data that's error encoded. 

    The encoding input is a linear array of length payloadSize * nCol bytes,
    where nRow is ideally set to the payload per strand. However, it's
    not a strict requirement.

    Bcol must be less than rs.field_charac, which is either 2**8 or 2**16.

    """
    def __init__(self,packetSize,errorSymbols,payloadSize,c_exp=8,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

        self._rs = get_reed_solomon(c_exp=c_exp)

        self._packetSize = packetSize
        self._errorSymbols = errorSymbols
        self._payloadSize = payloadSize
        
        self._lengthMessage = packetSize / payloadSize + errorSymbols

        assert(self._lengthMessage <= self._rs.field_charac) # requirement of RS
        
    """
    packet is a tuple with the key at position 0 and value at position 1. This should return
    the Reed-Solomon encoded byte array.
    """
    def _encode(self,packet):

        index = packet[0]
        data = packet[1][:]
        
        assert packet[0] < 256**3
        assert len(packet[1])<=self._packetSize

        # hopefully the last packet if it's not a full packetSize
        # so, pad it out with zeros to make an even number of strands if
        # isn't already
        
        if len(data) < self._packetSize:
            stats.inc("RSOuterCodec.padPacket")
            # normalize to multiple of payloadSize
            rem = len(data) % self._payloadSize
            data += [0]*(self._payloadSize-rem)

        #print data
        assert len(data) % self._payloadSize == 0
        rows = len(data) / self._payloadSize
        # construct message

        ecc = []        
        for i in range(self._payloadSize):
            message = data[i:len(data):self._payloadSize]
            mesecc = self._rs.rs_encode_msg(message, self._errorSymbols)
            ecc += mesecc[len(message):]
            #print "encode message=",mesecc
            #print message, ecc
            #print  len(message), len(ecc), self._errorSymbols
            #assert len(ecc) == self._errorSymbols

        tranpose = []
        for i in range( self._errorSymbols ):
            syms = ecc[i:len(ecc):self._errorSymbols]
            data += syms
                    
        # separate out key and value
        # convert mesecc into a string
        return packet[0],data

    """
    This function expects a list of unsigned integers in GF(256). For now, erasures are
    denoted with -1.
    """
    def _decode(self,packet):
        #print packet[1]
        data = [x for x in packet[1]]
        #print data
        rows = len(data) / self._payloadSize
        for i in range(self._payloadSize):
            message = data[i:len(data):self._payloadSize]            
            erasures = [ j for j in range(len(message)) if message[j]==-1 ]
            try:
                # correct the message
                corrected_message, corrected_ecc = \
                    self._rs.rs_correct_msg(message, \
                                            self._errorSymbols, \
                                            erase_pos=erasures)
                #print corrected_message
                #print corrected_ecc
                data[i:len(data):self._payloadSize] = corrected_message + corrected_ecc
                stats.inc("RSOuterCodec::correct")
            except ReedSolomonError as e:
                #print "couldn't correct block {}".format(message)
                stats.inc("RSOuterCodec::ReedSolomonError")
                # wrap exception into a library specific one when checking policy:
                wr_e = DNAReedSolomonOuterCodeError(msg=\
                             "RSOuterCodec found error at index={}".format(packet[0]))
                if self._Policy.allow(wr_e):
                    # fix -1
                    corrected_message = [ max(0,_) for _ in message ]
                    data[i:len(data):self._payloadSize] = corrected_message
                    pass # nothing to do
                else:
                    # policy doesn't allow handling, raise exception
                    raise wr_e
                # just proceed without further error correction
                pass
            except Exception as e:
                if self._Policy.allow(e):
                    pass
                else:
                    raise e

        # discard outer correction codes
        data = data[0:len(data)-self._errorSymbols*self._payloadSize]        
        return packet[0],data


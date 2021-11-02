#!/usr/bin/python
import os
import sys
import logging
logger = logging.getLogger("dna.storage.util.packetizedfile")
logger.addHandler(logging.NullHandler())





#converts a byte array into an integer, first byte is considered least significant 
def convertBytestoInt(l):
    _sum=0
    for i, val in enumerate(l):
        _sum+=ord(val)*(256**i)
    return _sum
"""
Write a file by receiving (index,value) packets of data of a specific size at a given index. Value
is of packetSize length.  index specifies where in the file to write the data.  Packets can be
received in any order. But, the file cannot be written until the file is complete, in other words,
it has received all size/packetSize packets.  The indices are assumed to be densely numbered from 0
to size/packetSize by 1.  Pythonically: range(0,size/packetSize+(size-size/packetSize*packetSize),1).
"""
class WritePacketizedFilestream:
    def __init__(self,fd,size,packetSize,minKey=0,zeroFillMissing=False):
        self.__fd = fd
        self.size = size
        self.__set_packetSize(packetSize)
        self.__data = {}
        self.minKey = minKey                      # inclusive
        self.zeroFillMissing = zeroFillMissing
        
    def has_key(self,key):
        return self.__data.has_key(key)
    def __setitem__(self,key,value):
        if (key >= self.minKey) and (key < self.maxKey):
            self.__data[key] = value
        #else:
            #print "not in range:",key,value,self.minKey,self.maxKey
    def __getitem__(self,key):
        assert key >= self.minKey and key < self.maxKey
        return self.__data[key]

    @property
    def maxKey(self):
        return self.minKey+self.numberOfPackets # exclusive
    
    @property
    def numberOfPackets(self):
        if (self.size % self.packetSize) > 0:
            return self.size / self.packetSize + 1
        else:
            return self.size / self.packetSize

    # packet property
    def __set_packetSize(self,val):
        self.__packetSize = val
    def __get_packetSize(self):
        return self.__packetSize
    packetSize = property(__get_packetSize,__set_packetSize)

    @property
    def lastPacketSize(self):
        # this function should never return 0, not the same as modulo
        # should return a number >= 1 and <= self.packetSize
        return self.size - (self.numberOfPackets-1)*self.packetSize

    @property
    def complete(self):
        if len(self.__data.keys()) < self.numberOfPackets:
            #print "too few keys {}".format(self.numberOfPackets)
            return False
        for i in range(self.minKey,self.maxKey):
            if not (i in self.__data.keys()):
                #print "missing key {}".format(i)
                return False
        return True

    def getMissingKeys(self):
        keys = self.__data.keys()
        missing = []
        for i in range(self.minKey,self.maxKey):
            if not i in keys:
                missing.append(i)
        return missing


    ## Warning: requires buffering the whole file!
    def write(self):
        emptyString = '\x00'*self.packetSize
        for i in range(self.minKey,self.maxKey):
            if self.__data.has_key(i):
                if i == self.maxKey-1:
                    self.__fd.write(self.__data[i][0:self.lastPacketSize])
                else:
                    self.__fd.write(self.__data[i])
            else:
                if self.zeroFillMissing:
                    self.__fd.write(emptyString)
                
        # for key,value in items:
        #     if i < key:
        #         while i < key:
        #             self.__fd.write(emptyString)
        #             i+=1
        #     if i == self.numberOfPackets-1:
        #         self.__fd.write(value[0:self.lastPacketSize])
        #     else:
        #         self.__fd.write(value)
        #     i+=1
        self.__fd.flush()
            
    def close(self):
        if self.__fd != sys.stdout and self.__fd != sys.stderr:
            self.__fd.close()

    def dummy_write(self):
        items = self.__data.items()
        items.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        i = 0
        emptyString = '\x00'*self.packetSize
        output_data=""
        for key,value in items:
            if value is not type(str):
                value2="".join(value)
            else:
                value2 = value
            if i < key:
                while i < key and i < self.numberOfPackets-1:
                    output_data+=emptyString
                    i+=1
            if i == self.numberOfPackets-1:
                output_data+=value2[0:self.lastPacketSize]
                i+=1
                break
            else:
                output_data+=value2
                i+=1
            if i==self.numberOfPackets:
                break
        while(i<self.numberOfPackets):
            if i == self.numberOfPackets-1:
                output_data+=emptyString[0:self.lastPacketSize]
            else:
                output_data+=emptyString
            i+=1
        #print i,self.numberOfPackets,len(output_data),self.size
        assert i == self.numberOfPackets and len(output_data)==self.size 
        return output_data



class WritePacketizedFile(WritePacketizedFilestream):
    def __init__(self,filename,size,packetSize):
        WritePacketizedFilestream.__init__(self,open(filename,"wb"),size,packetSize)
        self.__filename = filename



"""
Read a file by breaking it up into packetSize pieces. If the last piece is smaller
than packetSize, it's padded with zeros. Clients are guaranteed all packets are uniform
size.

This is an iterable object. Packets can be read using the iterator or by requesting
a specific index between 0 and size/packetSize.

This class expects a file descriptor. There is a derived class that accepts a filename.
"""
class ReadPacketizedFilestream:
    def __init__(self,fd):
        self.__fd = fd
        self.__set_packetSize(120) #packetSize is interpreted in terms of codewords
        self.__set_codewordSize(8) #codewordSize is interpreted in terms of bits, default to 8 bits so code that assumes byte arrays receive that
        assert self.codewordSize<32 #set a limit on how large a codeword can be 
        self.__read_size = 0
        self._RS=False
        self._bit_pointer=[0,0] #this pointer is going to be used to keep track of percisely where we are in the file when doing codeword sizes that are not multiples of bytes, [0] = byte [1] = bit offset in byte
    def read(self):
        
        #print "packet size {} byte pointer {} bit poiner {}".format(self.packetSize,self._bit_pointer[0],self._bit_pointer[1]) 
        bits_from_pointer=self.packetSize*self.codewordSize #total number of bits we should read
        #calculate the number of bytes to read 
        if ((self._bit_pointer[1]+bits_from_pointer) % 8) == 0:
            bytes_to_get=(self._bit_pointer[1]+bits_from_pointer)/8
        else:
            bytes_to_get=((self._bit_pointer[1]+bits_from_pointer)/8)+1
        
        #array of bytes to create codewords
        b = bytes(self.__fd.read(bytes_to_get))
        if not b: #done reading data from the file
            return b
        
        codeword_array=[]
        #extract values from the byte array into the codeword array, if we run out of bits, that is,8 bits is not a multiple of the codeword size, just pad it with 0's.
        bit_iter=self._bit_pointer[1]
        byte_iter=0
        for codeword_index in range(0,self.packetSize):
            #print "byte iter: {} bit_iter {}".format(byte_iter,bit_iter)
            assert byte_iter<len(b)
            if self.codewordSize + bit_iter <= 8: #bits come from the same byte
                mask=255>>(8-self.codewordSize)
                cw=(convertBytestoInt(b[byte_iter])>>(8-(bit_iter+self.codewordSize)))&mask 
                byte_iter=byte_iter+(bit_iter+self.codewordSize)/8
                bit_iter=(bit_iter+self.codewordSize)%8
                assert cw < 2**self.codewordSize
                #print "cw size {} cw {}".format(self.codewordSize,cw)
                codeword_array.append(cw)
            else: #bits come from different bytes
                built_cw=0x0000000000000000
                #need to know how many bytes the codeword spans
                if((self.codewordSize+bit_iter)%8==0):
                    bytes_touched=(self.codewordSize+bit_iter)/8
                    bits_to_last=8 #number of bits that come from the last byte
                else:
                    bytes_touched=((self.codewordSize+bit_iter)/8)+1
                    bits_to_last=((self.codewordSize+bit_iter)%8)
                bits_to_first=8-bit_iter #number of bits that come from the first byte
                #print "bit iter {}  bits to first {} bits to last {}  codeword size {}".format(bit_iter,bits_to_first,bits_to_last,self.codewordSize)
                assert (bits_to_first+bits_to_last+8*(bytes_touched-2)) == self.codewordSize #make sure that bits contributed from the first,last, and in between bytes, add to the codeword size
                for _ in range(0,bytes_touched):
                    #get the bits from each byte
                    if _ == 0:
                        mask=255>>(8-bits_to_first)
                        #print "first b[byte_iter+_] {}".format(bin(convertBytestoInt(b[byte_iter+_])))
                        partial_cw=(convertBytestoInt(b[byte_iter+_])>>((8-(bit_iter+bits_to_first))))&mask
                        difference=self.codewordSize-bits_to_first
                        built_cw=(partial_cw)<<difference
                        #print difference
                        #print"partial  cw {} first cw after {}".format(bin(partial_cw),bin(built_cw))
                    elif _ == bytes_touched-1:
                        mask=255>>(8-bits_to_last)
                        #print "last b[byte_iter+_] {}".format(bin(convertBytestoInt(b[byte_iter+_])))                                        
                        partial_cw=(convertBytestoInt(b[byte_iter+_])>>((8-bits_to_last)))&mask
                        built_cw=built_cw|partial_cw
                        #print"last>> partial cw {} last cw after {}".format(bin(partial_cw),bin(built_cw))
                    else:
                        partial_cw=convertBytestoInt(b[byte_iter+_])
                        built_cw=built_cw|(partial_cw<<(self.codewordSize-(bits_to_first+_*8)))
                        #print"partial cw {} middle cw after {}".format(bin(partial_cw),bin(built_cw))
                    if byte_iter+(_+1)>=len(b): break
                assert built_cw < 2**self.codewordSize
                codeword_array.append(built_cw)
                byte_iter+=(self.codewordSize+bit_iter)/8
                bit_iter=(self.codewordSize+bit_iter)%8
            if byte_iter>=len(b): break #past the number of bytes we have 
                
        self.__read_size += len(codeword_array)*self.codewordSize #read size will be the number of bits that we extract out of the byte array
        
        #only pad out if it is not RS, RS encodes by blocks so padding may lead to many useless strands
        if b  and len(codeword_array) != self.packetSize and not self._RS:
            codeword_array+=[1 for _ in range(0,(self.packetSize-len(codeword_array)))]
            assert len(codeword_array)==self.packetSize

        self._bit_pointer[1]=(self._bit_pointer[1]+bits_from_pointer)%8 #increase bit pointer
        if self._bit_pointer[1]%8==0:
            self._bit_pointer[0]+=bytes_to_get # increase byte pointer by the number of bytes we got since the starting bit position + #bits we got is on a byte boundary
        else:
            self._bit_pointer[0]+=(bytes_to_get-1)#dont want to seek past a byte we have not finished
        self.__fd.seek(self._bit_pointer[0])#reset the byte pointer for the next read(), makes sure we read the bits that may be let over from a byte
        
        return codeword_array

    # packet property
    def __set_packetSize(self,val):
        self.__packetSize = val
    def __get_packetSize(self):
        return self.__packetSize
    def __set_codewordSize(self,val):
        self.__codewordSize=val
    def __get_codewordSize(self):
        return self.__codewordSize
    
    packetSize = property(__get_packetSize,__set_packetSize)
    codewordSize=property(__get_codewordSize,__set_codewordSize)

    # iterator
    def __iter__(self):
        # reset to begnning of file
        self.__fd.seek(0)
        self._bit_pointer=[0,0]
        return self

    def next(self):
        b = self.read()
        if b:
            return b
        else:
            raise StopIteration()

    def __getitem__(self,key):
        self.__fd.seek(key*self.packetSize)
        return self.read()

    @property
    def size(self):
        return os.fstat(self.__fd.fileno()).st_size

    @property
    def bytes_read(self):
        if self.__read_size % 8:
            return self.__read_size/8
        else:
            return self.__read_size/8+1

    @property
    def bits_read(self):
        return self.__read_size
    
    @property
    def numberOfPackets(self):
        if ((self.size*8) % (self.packetSize*self.codewordSize))==0:
            return (self.size*8) / (self.packetSize*self.codewordSize) #convert to bits for both to normalize
        else:
            return ((self.size*8) / (self.packetSize*self.codewordSize)) + 1


class ReadPacketizedFile(ReadPacketizedFilestream):
    def __init__(self,filename):
        ReadPacketizedFilestream.__init__(self,open(filename,"rb"))
        self.__filename = filename
    @property
    def filename(self):
        return self.__filename


if __name__ == "__main__":
    import os
    import sys
    from random import randint
    packetFile = ReadPacketizedFile(sys.argv[1])

    out = WritePacketizedFile("output.d",packetFile.size,120)

    assert out.complete==False

    i = 0
    for p in packetFile:
        out[i] = p
        i += 1

    assert out.complete==True
    out.write()

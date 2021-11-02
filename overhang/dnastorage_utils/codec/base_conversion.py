 #!/usr/bin/python
from copy import copy
from random import *

bases = ['A', 'C', 'G', 'T']

values = { 'A' : 0 ,
           'C' : 1 ,
           'G' : 2 ,
           'T' : 3   }

def randomTernary(length):
    return "".join([ bases[x%3] for x in range(length) ])

def convertTernaryHelper(dec,s):
    m = dec % 3
    q = dec / 3
    s = s + bases[m]
    if q > 0:
        return convertTernaryHelper(q,s)
    else:
        return s

def convertTernary(dec,length):
    s = convertTernaryHelper(dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertQuarnaryHelper(dec,s):
    m = dec % 4
    q = dec / 4
    s = s + bases[m]
    if q > 0:
        return convertQuarnaryHelper(q,s)
    else:
        return s

def convertQuarnary(dec,length):
    s = convertQuarnaryHelper(dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertBaseHelper(base,dec,s):
    m = dec % base
    q = dec / base
    #print(s)
    s = s + bases[m]
    #print(s)
    if q > 0:
        return convertBaseHelper(base,q,s)
    else:
        return s
    
def convertBytetoBinary(x,length,reverse=False):
    s=convertBytetoBinaryHelper(x,'')
    s=s.ljust(length,"0")
    #print s
    if reverse == True:
        return s[::-1]
    else:
        return s

def convertBytetoBinaryHelper(x,s): #convert byte x to a binary string
    binary=['0','1']
    m=x%2
    q=x/2
    s=s+binary[m]
    if q >0:
        return convertBytetoBinaryHelper(q,s)
    else:
        return s

def convertBase(base,dec,length):
    s = convertBaseHelper(base,dec,'')
    s = s.ljust(length,bases[0])
    return s

def convertFromBase(base,s):
    val = 0
    power = 1
    for i in range(len(s)):
        val += values[s[i]]*power
        power *= base
    return val

def convertIntToBytes(val,num_bytes):
    val = int(val)
    if val == None:
        return [-1 for _ in range(num_bytes)]
    else:
        #print ("HERE!",val,num_bytes,range(num_bytes))
        l = [(val & (0xff << pos*8)) >> pos*8 for pos in range(num_bytes)]
        return l

def convertIntToCodeword(val,num_codewords,codewordSize=8): #this function takes an integer and creates a string of codewords 
    val=int(val)
    if val==None:
        return [-1 for _ in range(num_codewords)]
    else:
        mask=0xffffffff>>(32-codewordSize) #make sure we only get codeword sized chunks, make sure the mask is large enough to support larger codewords 
        l=[(val&(mask<<pos*codewordSize))>>pos*codewordSize for pos in range(num_codewords)]
        return l

    
def convertBytesToInt(l):
    _sum = 0
    for i,val in enumerate(l):
        _sum += val * (256**i)
    return _sum

ibases = ['A', 'C', 'T']

ivalues = { 'A' : 0 ,
           'C' : 1 ,
           'G' : 3 ,
           'T' : 2   }

def encodeWithExclusionHelper(dec,s,excluded):
    if len(s) < len(excluded):
        b = len(excluded[len(s)])
        m = dec % b
        q = dec / b
        s += excluded[len(s)][m]
    else:
        m = dec % 3
        q = dec / 3
        s += ibases[m]

    if q > 0 or len(s) < len(excluded):
        return encodeWithExclusionHelper(q,s,excluded)
    else:
        return s

def encodeWithExclusion(dec,length,primer):
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i] != 'G':
            b.remove(plist[i])
        excluded.append(b)
    #print excluded
    s = encodeWithExclusionHelper(dec,"",excluded)
    if len(s) < len(primer):
        for i in range(len(s),min(length,len(primer))):
            s += excluded[i][0]
    s = s.ljust(length,ibases[0])
    return s


def decodeWithExclusion(s,primer):
    val = 0
    power = 1
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i]!='G':
            b.remove(plist[i])
        excluded.append( { b[i] : i for i in range(len(b)) }  )
    #print excluded
    base = 3
    for i in range(len(s)):
        if i < len(excluded):
            base = len(excluded[i])
            val += excluded[i][s[i]] * power
        else:
            base = 3
            val += ivalues[s[i]] * power
        power *= base
    return val


if __name__ == "__main__":
    import math

    primer = "GTCTCGTGGGCTCGG"
    #for i in range(3**5):
    #    s = encodeWithExclusion(i,10,primer)
    #    print "{} = {}  ({})".format(i,s,decodeWithExclusion(s,primer))
    #    assert i == decodeWithExclusion(s,primer)

    for i in range(2**8):
        print (convertBase(2,i,8))

    for i in range (len(primer),60):
        val = (2**len(primer)) * (3 **(i-len(primer)))
        print ("Max({}) = {}".format(i,math.log(val,2)))

    for i in [2**64, 2**74, 2**80]:
        s = encodeWithExclusion(i,100,primer)
        print (s)

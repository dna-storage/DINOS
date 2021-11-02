from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.arch.strand import *
from dnastorage.arch.builder import customize_RS_CFC8
from dnastorage.arch.builder import build_overhang_bitstring_strand
from dnastorage.exceptions import *


def ENC_OH_BITSTRING_XXX(pf, primer5, primer3,strand_length=512,num_overhangs=3,bits_per_block=1): #encoding used to experimentally evaluate overhang construction
    #strand length is the length of a strand in codewords, and bits_per_block is the size of the codeword in bits
    enc = build_overhang_bitstring_strand(True,pf,primer5,primer3,strand_length,num_overhangs,bits_per_block)
    return enc

def DEC_OH_BITSTRING_XXX(pf, primer5,primer3,strand_length,num_overhangs,bits_per_block):
    assert 0 #not implemented yet

FileSystemFormats = {
    0x3000 : [0x3000, "Variable", "Variable", "OH_BITSTRING_XXX", "Experimental encoding to evaluate overhang construction",ENC_OH_BITSTRING_XXX,DEC_OH_BITSTRING_XXX],
}


def file_system_formats():
    return [ v[3] for k,v in FileSystemFormats.items() ]

_abbrevFileSystemDict = { v[3] : v for k,v in FileSystemFormats.items() }

def file_system_format_description(formatid):
    return FileSystemFormats[formatid][4]

def file_system_format_packetsize(formatid):
    return FileSystemFormats[formatid][2]

def file_system_encoder(formatid):
    return FileSystemFormats[formatid][5]

def file_system_decoder(formatid):
    return FileSystemFormats[formatid][6]

def file_system_encoder_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][5]

def file_system_decoder_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][6]

def file_system_formatid_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][0]




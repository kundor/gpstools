#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:58:56 2017

@author: nima9589
"""
from functools import reduce
from operator import xor
from hashlib import md5
from binascii import crc32, crc_hqx
import io

SYNC = b'\xC2'

RECS = {192 : 'SNR',
        193 : 'Housekeeping',
        194 : 'Meta'}

def info(*args, logfile=open('/tmp/loggin', 'a')):
    print(*args, file=logfile, flush=True)
    

def readsync(strm):
    """Read a sync byte from the stream. If we don't see one, scan forward for it."""
    read = 0
    while strm.read(1) != SYNC:
        read += 1
    if read:
        info("Skipped", str(read), "bytes before SYNC.")


def read_ubnxi(strm):
    """Read a little-endian ubnxi (unsigned binex integer) from the stream.

    Consume one to four bytes from the stream, stopping at any with the high bit set.
    High bits in any but the last byte are ignored. The constituent bytes and
    converted int are returned.
    """
    val = 0
    byts = b''
    for b in range(3):
        byts += strm.read(1)
        byt = byts[-1] # indexing bytes converts to int
        if b < 3:
            masked = byt & 127
        else:
            masked = byt
        val += masked << (7*b)
        if byt < 128:
            break
    if b and val < 2**(7*b):
        raise ValueError('Overlong ubnxi encoding')
    return byts, val

def checksum(msg):
    """Compute BINEX checksum (regular CRC) of the given bytes."""
    if len(msg) < 128:
        # 1-byte checksum: 8-bit XOR of all bytes
        return reduce(xor, msg).to_bytes(1, 'little')
    elif len(msg) < 4096:
        # 2-byte CRC, generating polynomial x^16 + x^12 + x^5 + x^0
        # AKA CRC-16-CCITT
        return crc_hqx(msg, 0).to_bytes(2, 'little')
    elif len(msg) < 1048576:
        # 4-byte CRC, generating polynomial x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x+1
        # this is the same as zlib's crc32 function
        return crc32(msg).to_bytes(4, 'little')
    else:
        # 16-byte CRC: 128-bit MD5 checksum
        return md5(msg).digest()

def peeker(strm, num):
    """Read num bytes from stream, without advancing."""
    vals = strm.peek(num)[:num]
# BufferedReader.peek just returns the entire buffer -- slice it down to size
    assert len(vals) == num, "FIXME: Need to read more into buffer to peek ahead"
    return vals

def verify(strm, msg):
    """Compute the checksum for msg and match it in the stream.
    
    The correct number of bytes for the length of msg are consumed from the stream
    (0--127 bytes: 1 byte read; 128--4095: 2 read; 4096--1048575: 4 read; otherwise 16).
    If they do not match, we just complain.
    """
    compcrc = checksum(msg)
    readcrc = strm.read(len(compcrc))
    assert len(readcrc) == len(compcrc)
    if compcrc != readcrc:
        info('Bad checksum! Computed', compcrc, '\n'
             '              Found   ', readcrc)

`
        

        
def read_record(strm):
    readsync(strm)
    idbytes, recid = read_ubnxi(strm)
    if recid not in RECS:
        info('Record ID', recid, 'unknown')
    lenbytes, msglen = read_ubnxi(strm)
    msg = strm.read(msglen)
    assert len(msg) == msglen
    verify(strm, idbytes + lenbytes + msg)
    

"""
Functions to process little-endian, forward-readable, regular-CRC BINEX files.
See http://binex.unavco.org/ for more information.
"""
import time
from functools import reduce
from operator import xor
from hashlib import md5
from binascii import crc32, crc_hqx
from utility import info, debug

SYNC = b'\xC2'

def readsync(strm):
    """Read a sync byte from the stream. If we don't see one, scan forward for it."""
    read = b''
    r = strm.read(1)
    while r != SYNC:
        if not r:
            raise EOFError('End of file encountered before sync byte')
        read += r
        r = strm.read(1)
    if read:
        info("Skipped", len(read), "bytes before SYNC:", read.hex())


def read_ubnxi(strm):
    """Read a little-endian ubnxi (unsigned binex integer) from the stream.

    Consume one to four bytes from the stream, stopping at any with the high bit set.
    High bits in any but the last byte are ignored. The constituent bytes and
    converted int are returned.
    """
    val = 0
    byts = b''
    for b in range(4):
        byt = strm.read(1)
        if not byt:
            raise EOFError('ubnxi cut off by end of stream')
        byts += byt
        byt = ord(byt)
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

def bnx_prn(byt, badprns):
    """Interpret a BINEX 1-byte satellite ID.

    See http://binex.unavco.org/conventions.html#SVid1_details
    Return PRN for GPS (1--32), SBAS (120--151), QZSS (193--200),
    32 + PRN for Galileo (33--64), and 64 + slot # for GLONASS (65--96).
    """
    if isinstance(byt, bytes):
        byt = ord(byt)
    satsys = byt >> 5
    if satsys:
        info("Non-GPS satellite system:", satsys)
        if satsys > 5:
            raise ValueError('   Impossible satsys! Greater than 5 reserved.')
    # FIXME: KLUDGE: currently GPS PRNs are not being subtraced by 1 at the source
    #  (first value below should be 1)
    corrections = [1 - badprns, 65, 120, 33, 201, 193, 0, 0]
    # GPS, add 1 for PRN
    # GLONASS, add 1 for slot #, return 64 + slot #
    # SBAS, add 120 for PRN
    # Galileo, add 1 for prn, and 32 to map into 33--64
    # Beidou / COMPASS: undefined, but NMEA reserves 201--235
    # QZSS, add 193 for PRN
    return (byt & 31) + corrections[satsys]

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

def binex_to_weeksow(msg):
    """Convert 6-byte (millisecond) binex timestamp to a pair: GPS week, second of week.

    See http://binex.unavco.org/#time_stamps.
    """
    gmin = int.from_bytes(msg[:4], 'little')
    gmsc = int.from_bytes(msg[4:6], 'little')
    return gmin//10080, (gmin % 10080)*60 + gmsc/1000

def snr_msg(msg):
    """Read SNR record. Return tuple of RXID, (GPS week, second of week) pair, and
    a list of PRN, SNR pairs."""
    rxid = msg[0]
    weeksow = binex_to_weeksow(msg[1:7])
    badprns = (rxid == 5 and weeksow[0] <= 1940 and weeksow[1] < 2*24*60*60) # Incorrect prns were put out in this software version
    snrs = [(bnx_prn(msg[i], badprns), int.from_bytes(msg[i+1:i+3], 'little'))
            for i in range(7, len(msg), 3)]
    return rxid, weeksow, snrs

def hk_msg(msg):
    """Read HK record. Return tuple of RXID, (GPS week, second of week) pair,
    MAC, longitude, latitude, altitude, voltage, temperature, message count, and
    error flags."""
    uint = lambda x: int.from_bytes(x, 'little')
    sint = lambda x: int.from_bytes(x, 'little', signed=True)
    fields = ((0, 1, uint),   # RXID
              (1, 7, binex_to_weeksow), # time
              (7, 15, uint),  # MAC
              (15, 19, sint), # lon, decimal degrees * 10**7
              (19, 23, sint), # lat, decimal degrees * 10**7
              (23, 27, sint), # alt, millimeters above ellipsoid
              (27, 29, uint), # volt
              (29, 30, uint), # temp, only top byte used for now. degrees celsius
              (31, 33, uint), # msgct
              (33, 34, uint)) # error flags
    return [fn(msg[a:b]) for a,b,fn in fields]

RECS = {192 : snr_msg,
        193 : hk_msg}

def read_record(strm):
    readsync(strm)
    idbytes, recid = read_ubnxi(strm)
    lenbytes, msglen = read_ubnxi(strm)
    if msglen > 128:
        # We have no messages so long. Probably this is corrupted.
        # Instead of reading further, raise a ValueError, so the enclosing
        # loop can try for the next sync byte
        raise ValueError('Message too long (id {}, length {})'.format(recid, msglen))
    # FIXME: KLUDGE: currently length is wrong in the stream for SNR records
    if recid == 192:
        corlen = round(msglen/3)*3 + 1 # should always be 1 mod 3
        if corlen != msglen:
            debug('Correcting SNR message length from', msglen, 'to', corlen)
            msglen = corlen
    msg = strm.read(msglen)
    tries = 0
    while len(msg) < msglen:
        tries += 1
        if tries > 5:
            break
        time.sleep(.1)
        msg += strm.read(msglen - len(msg))
    if len(msg) != msglen:
        raise ValueError('Could not read full message (received {} bytes out of {})'.format(len(msg), msglen))
    verify(strm, idbytes + lenbytes + msg)
    # When verification fails, we could try advancing to the next sync byte
    # within the already-read message, and attempt to read a record there.
    # (Currently, we start trying after the end of all data previously read.)
    if recid not in RECS:
        info('Record ID', recid, 'unknown')
        return recid, msg
    return recid, RECS[recid](msg)

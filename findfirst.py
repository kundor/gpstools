# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 09:20:28 2016

@author: nima9589
"""
from numba import jit

@jit(nopython=True)
def findfirstgt(vec, val):
    """Find index of first value in vec greater than val, or -1 if none exist."""
    for i in range(len(vec)):
        if vec[i] > val:
            return i
    return -1

@jit(nopython=True)
def findfirstgte(vec, val):
    """Index of first value in vec greater than or equal to val, or -1 if none exist."""
    for i in range(len(vec)):
        if vec[i] >= val:
            return i
    return -1

@jit(nopython=True)
def findfirstlt(vec, val):
    """Find index of first value in vec less than val, or -1 if none exist."""
    for i in range(len(vec)):
        if vec[i] < val:
            return i
    return -1

@jit(nopython=True)
def findfirstlte(vec, val):
    """Index of first value in vec less than or equal to val, or -1 if none exist."""
    for i in range(len(vec)):
        if vec[i] <= val:
            return i
    return -1

def makefindpred(pred):
    """Create a JIT function to find first value satisfying a given predicate.

    The returned function takes an array and returns the index of the first value
    satisfying pred, which must itself be a JIT function.
    """
# This workaround is necessary because nopython JIT functions cannot take functions
# as arguments.
    @jit(nopython=True)
    def findfirstpred(vec):
        for i in range(len(vec)):
            if pred(vec[i]):
                return i
        return -1
    return findfirstpred


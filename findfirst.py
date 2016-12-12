#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 09:20:28 2016

@author: nima9589
"""
from numba import jit

@jit(nopython=True)
def findfirstgt(vec, val):
    for i in range(len(vec)):
        if vec[i] > val:
            return i
    return -1

@jit(nopython=True)
def findfirstgte(vec, val):
    for i in range(len(vec)):
        if vec[i] >= val:
            return i
    return -1

@jit(nopython=True)
def findfirstlt(vec, val):
    for i in range(len(vec)):
        if vec[i] < val:
            return i
    return -1

@jit(nopython=True)
def findfirstlte(vec, val):
    for i in range(len(vec)):
        if vec[i] <= val:
            return i
    return -1

def makefindpred(pred):
    @jit(nopython=True)
    def findfirstpred(vec):
        for i in range(len(vec)):
            if pred(vec[i]):
                return i
        return -1
    return findfirstpred


#!/usr/bin/python
__author__ = 'Josh'
import cPickle

def make_pickle(obj, filename):
    with open(filename, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
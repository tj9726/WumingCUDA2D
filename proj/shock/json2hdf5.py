#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" JSON to HDF5 converter
"""

import os
import sys
import numpy as np
import json
import h5py

_DEBUG = True


def json2hdf5(jsonfile, hdffile=None, verbose=True):
    if hdffile is None:
        hdffile = os.path.splitext(jsonfile)[0] + '.h5'

    if verbose:
        print('Processing {} to produce {}'.format(jsonfile, hdffile))

    # read json and store in memory
    with open(jsonfile, 'r') as fp:
        obj = json.load(fp)

    # meta data
    meta = obj.get('meta')

    # check endian
    endian = meta.get('endian')
    if endian == 1:          # little endian
        byteorder = '<'
    elif endian == 16777216: # big endian
        byteorder = '>'
    else:
        errmsg = 'unrecognized endian flag: {}'.format(endian)
        report_error(errmsg)

    # check raw data file
    rawfile = meta.get('rawfile')
    if not os.path.exists(rawfile):
        errmsg = 'rawfile {} does not exist'.format(rawfile)
        report_error(errmsg)

    #
    # create hdf5
    #
    with h5py.File(hdffile, 'w') as h5fp, \
         open(rawfile, 'r') as rawfp:
        # attribute
        attribute = obj.get('attribute', [])
        for name in attribute:
            attr = attribute.get(name)
            data = read_data(rawfp, attr, byteorder)
            h5fp.attrs.create(name, data)
            if verbose:
                print('  - attribute "{}" has been created '.format(name), end='')
                print('with data : {}'.format(data))

        # dataset
        dataset = obj.get('dataset', [])
        for name in dataset:
            data = dataset.get(name)
            offset, dsize, dtype, shape = read_info(data, byteorder)
            shape = shape[::-1] # assume column-major shape
            ext = ((rawfile, offset, dsize), )
            h5fp.create_dataset(name, shape=shape, dtype=dtype, external=ext)
            if verbose:
                print('  - dataset "{}" has been created '.format(name), end='')
                print('with dtype = "{}" and shape = "{}"'.format(dtype, shape))

    if verbose:
        print('done !')


def read_data(fp, obj, byteorder):
    offset   = obj['offset']
    datatype = byteorder + obj['datatype']
    shape    = obj['shape']
    size     = np.product(shape)
    fp.seek(offset)
    x = np.fromfile(fp, datatype, size).reshape(shape)
    if len(shape) == 1 and shape[0] == 1:
        x = x[0]
    return x


def read_info(obj, byteorder):
    offset   = obj['offset']
    datatype = byteorder + obj['datatype']
    shape    = obj['shape']
    datasize = np.product(shape) * np.dtype(datatype).itemsize
    return offset, datasize, datatype, shape


def error_report(msg):
    print('Error: {}'.format(msg))
    if _DEBUG:
        raise ValueError(errmsg)


if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('-v', '--verbose', dest='verbose',
                      action='store_true', default=False,
                      help='print verbose messages')
    options, args = parser.parse_args()
    for f in args:
        if os.path.exists(f) and os.path.splitext(f)[1] == '.json':
            json2hdf5(f, verbose=options.verbose)

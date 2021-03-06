# -*- coding: utf-8 -*-
# Copyright 2014 by Forschungszentrum Juelich GmbH
# Author: J. Caron
#
"""This module provides the :class:`~.SemperFormat` class which allows reading and writing `.unf`
files which are used in Semper."""


import os
from time import strftime
import numpy as np
import struct

import logging


__all__ = ['SemperFormat']


class SemperFormat(object):

    '''Class for importing and exporting Semper `.unf`-files.

    The :class:`~.SemperFormat` class represents a Semper binary file format with a header, which
    holds additional information. Floating point images can be converted to :class:`~.PhaseMaps',
    complex value images can be converted to :class:`~.MagData'-slices (2D magnetization along
    a simple axis) and vice versa with the corresponding methods and classmethods. `.unf`-files
    can also be saved and read from files.

    Attributes
    ----------
    dim : tuple (N=3)
        Dimensions of the data.
    a : float
        Grid spacing (nm per pixel).
    data : :class:`~numpy.ndarray` (N=3)
        The phase map or magnetization information in a 3D array (with one slice).
    iclass : int
        Defines the image class defined in `ICLASS_DICT`. Normally `image` (1) is chosen.
    iform : int
        Defines the data format defined in 'IFORM_DICT'.
    iversn :
        Current `.unf`-format version. Current: 2.
    ilabel : int
        Defines if a label is present (1) or not (0).
    iformat : int
        Defines if the file is formatted (1) or not (0).
    title : string
        Title of the file (not to be confused with the filename).
    iwp : int
        Write protect flag, determining if picture is (1) or is not (0) write-projtected.
    ipltyp : int
        Position type list. Standard seems to be 0 (picture not a position list).
    date : string
        The date of file construction.
    iccoln : int
        Column number of picture origin.
    icrown : int
        Row number of picture origin.
    iclayn : int
        Layer number of picture origin.

    ''' # TODO: Docstrings (everywhere)!

    _log = logging.getLogger(__name__)

    ICLASS_DICT = {1: 'image', 2: 'macro', 3: 'fourier', 4: 'spectrum',
                   5: 'correlation', 6: 'undefined', 7: 'walsh', 8: 'position list',
                   9: 'histogram', 10: 'display look-up table'}

    ICLASS_DICT_INV = {v: k for k, v in ICLASS_DICT.iteritems()}

    IFORM_DICT = {0: np.byte, 1: np.int32, 2: np.float32, 3: np.complex64}

    IFORM_DICT_INV = {v: k for k, v in IFORM_DICT.iteritems()}

    def __init__(self, arg_dict):
        self._log.debug('Calling __init__')
        self.data = arg_dict['data']
        self.title = arg_dict['title']
        self.offsets = arg_dict['offsets']
        self.scales = arg_dict['scales']
        self.units = arg_dict['units']
        self.date = arg_dict['date']
        self.iclass = arg_dict['ICLASS']
        self.iform = arg_dict['IFORM']
        self.iversn = arg_dict['IVERSN']
        self.ilabel = arg_dict['ILABEL']
        self.iformat = arg_dict['IFORMAT']
        self.iwp = arg_dict['IWP']
        self.ipltyp = arg_dict['IPLTYP']
        self.iccoln = arg_dict['ICCOLN']
        self.icrown = arg_dict['ICROWN']
        self.iclayn = arg_dict['ICLAYN']
        self._log.debug('Created '+str(self))

    @classmethod
    def from_file(self, filename):
        '''Load a `.unf`-file into a :class:`~.SemperFormat` object.

        Parameters
        ----------
        filename : string
            The name of the unf-file from which to load the data. Standard format is '\*.unf'.

        Returns
        -------
        semper : :class:`~.SemperFormat` (N=1)
            Semper file format object containing the loaded information.

        '''
        self._log.debug('Calling from_file')
        # Construct path if filename isn't already absolute:
        if not os.path.isabs(filename):
            from pyramid import DIR_FILES
            directory = os.path.join(DIR_FILES, 'semper')
            if not os.path.exists(directory):
                os.makedirs(directory)
            filename = os.path.join(directory, filename)
        with open(filename, 'rb') as f:
            # Read header:
            rec_length = np.frombuffer(f.read(4), dtype=np.int32)[0]  # length of header
            header = np.frombuffer(f.read(rec_length), dtype=np.int16)
            ncol, nrow, nlay = header[:3]
            iclass = header[3]
            iform = header[4]
            data_format = self.IFORM_DICT[iform]
            iflag = header[5]
            iversn, remain = divmod(iflag, 10000)
            ilabel, ntitle = divmod(remain, 1000)
            iformat = header[6] if len(header) == 7 else None
            assert np.frombuffer(f.read(4), dtype=np.int32)[0] == rec_length
            # Read title:
            title = ''
            if ntitle > 0:
                assert np.frombuffer(f.read(4), dtype=np.int32)[0] == ntitle  # length of title
                title_bytes = np.frombuffer(f.read(ntitle), dtype=np.byte)
                title = ''.join(map(chr, title_bytes))
                assert np.frombuffer(f.read(4), dtype=np.int32)[0] == ntitle
            # Read label:
            iwp, date, range_string, ipltype, a = [None] * 5  # Initialization!
            iccoln, icrown, iclayn = [None] * 3
            if ilabel:
                rec_length = np.frombuffer(f.read(4), dtype=np.int32)[0]  # length of label
                label = np.frombuffer(f.read(512), dtype=np.int16)
                assert ''.join([chr(l) for l in label[:6]]) == 'Semper'
                assert struct.unpack('>h', ''.join([chr(x) for x in label[6:8]]))[0] == ncol
                assert struct.unpack('>h', ''.join([chr(x) for x in label[8:10]]))[0] == nrow
                assert struct.unpack('>h', ''.join([chr(x) for x in label[10:12]]))[0] == nlay
                iccoln = struct.unpack('>h', ''.join([chr(x) for x in label[12:14]]))[0]
                icrown = struct.unpack('>h', ''.join([chr(x) for x in label[14:16]]))[0]
                iclayn = struct.unpack('>h', ''.join([chr(x) for x in label[16:18]]))[0]
                assert label[18] == iclass
                assert label[19] == iform
                iwp = label[20]
                date = '{}-{}-{} {}:{}:{}'.format(label[21]+1900, *label[22:27])
                # No test for ncrang, range is extracted from data itself (also prone to errors)!
                ipltyp = label[55]  # position list type
                real_coords = label[62]
                dz, dy, dx, z0, y0, x0 = [1., 1., 1., 0., 0., 0.]
                if real_coords:
                    dz = struct.unpack('<f', ''.join([chr(x) for x in label[75:79]]))[0]
                    z0 = struct.unpack('<f', ''.join([chr(x) for x in label[79:83]]))[0]
                    dy = struct.unpack('<f', ''.join([chr(x) for x in label[83:87]]))[0]
                    y0 = struct.unpack('<f', ''.join([chr(x) for x in label[87:91]]))[0]
                    dx = struct.unpack('<f', ''.join([chr(x) for x in label[91:95]]))[0]
                    x0 = struct.unpack('<f', ''.join([chr(x) for x in label[95:99]]))[0]
                assert ''.join([str(unichr(l)) for l in label[100:100+ntitle]]) == title
                ux = ''.join([chr(l) for l in label[244:248]]).replace('\x00', '')
                uy = ''.join([chr(l) for l in label[248:252]]).replace('\x00', '')
                uz = ''.join([chr(l) for l in label[252:256]]).replace('\x00', '')
                assert np.frombuffer(f.read(4), dtype=np.int32)[0] == rec_length
            # Read picture data:
            data = np.empty((nlay, nrow, ncol), dtype=data_format)
            for k in range(nlay):
                for j in range(nrow):
                    rec_length = np.frombuffer(f.read(4), dtype=np.int32)[0]  # length of row
                    row = np.frombuffer(f.read(rec_length), dtype=data_format)
                    data[k, j, :] = row
                    assert np.frombuffer(f.read(4), dtype=np.int32)[0] == rec_length
        arg_dict = {}
        arg_dict['data'] = data
        arg_dict['title'] = title
        arg_dict['offsets'] = (x0, y0, z0)
        arg_dict['scales'] = (dx, dy, dz)
        arg_dict['units'] = (ux, uy, uz)
        arg_dict['date'] = date
        arg_dict['ICLASS'] = iclass
        arg_dict['IFORM'] = iform
        arg_dict['IVERSN'] = iversn
        arg_dict['ILABEL'] = ilabel
        arg_dict['IFORMAT'] = iformat
        arg_dict['IWP'] = iwp
        arg_dict['IPLTYP'] = ipltyp
        arg_dict['ICCOLN'] = iccoln
        arg_dict['ICROWN'] = icrown
        arg_dict['ICLAYN'] = iclayn
        return SemperFormat(arg_dict)

    def to_file(self, filename='semper.unf', skip_header=False):
        '''Save a :class:`~.SemperFormat` to a file.

        Parameters
        ----------
        filename : string, optional
            The name of the unf-file to which the data should be written.
        skip_header : boolean, optional
            Determines if the header, title and label should be skipped (useful for some other
            programs). Default is False.

        Returns
        -------
        None

        '''
        self._log.debug('Calling to_file')
        # Construct path if filename isn't already absolute:
        if not os.path.isabs(filename):
            from pyramid import DIR_FILES
            directory = os.path.join(DIR_FILES, 'semper')
            if not os.path.exists(directory):
                os.makedirs(directory)
            filename = os.path.join(directory, filename)
        nlay, nrow, ncol = self.data.shape
        with open(filename, 'wb') as f:
            if not skip_header:
                # Create header:
                header = []
                header.extend(reversed(list(self.data.shape)))  # inverse order!
                header.append(self.iclass)
                header.append(self.iform)
                header.append(self.iversn*10000 + self.ilabel*1000 + len(self.title))
                if self.iformat is not None:
                    header.append(self.iformat)
                # Write header:
                f.write(struct.pack('I', 2*len(header)))  # record length, 4 byte format!
                for element in header:
                    f.write(struct.pack('h', element))  # 2 byte format!
                f.write(struct.pack('I', 2*len(header)))  # record length!
                # Write title:
                f.write(struct.pack('I', len(self.title)))  # record length, 4 byte format!
                f.write(self.title)
                f.write(struct.pack('I', len(self.title)))  # record length, 4 byte format!
                # Create label:
                if self.ilabel:
                    label = np.zeros(256, dtype=np.int32)
                    label[:6] = [ord(c) for c in 'Semper']
                    label[6:8] = divmod(ncol, 256)
                    label[8:10] = divmod(nrow, 256)
                    label[10:12] = divmod(nlay, 256)
                    label[12:14] = divmod(self.iccoln, 256)
                    label[14:16] = divmod(self.icrown, 256)
                    label[16:18] = divmod(self.iclayn, 256)
                    label[18] = self.iclass
                    label[19] = self.iform
                    label[20] = self.iwp
                    year, time = self.date.split(' ')
                    label[21:24] = map(int, year.split('-'))
                    label[21] -= 1900
                    label[24:27] = map(int, time.split(':'))
                    range_string = '{:.4g},{:.4g}'.format(self.data.min(), self.data.max())
                    ncrang = len(range_string)
                    label[27] = ncrang
                    label[28:28+ncrang] = [ord(c) for c in range_string]
                    label[55] = self.ipltyp
                    label[62] = 1  # Use real coords!
                    label[75:79] = [ord(c) for c in struct.pack('<f', self.scales[2])]  # DZ
                    label[79:83] = [ord(c) for c in struct.pack('<f', self.offsets[2])]  # Z0
                    label[83:87] = [ord(c) for c in struct.pack('<f', self.scales[1])]  # DY
                    label[87:91] = [ord(c) for c in struct.pack('<f', self.offsets[1])]  # Y0
                    label[91:95] = [ord(c) for c in struct.pack('<f', self.scales[0])]  # DX
                    label[95:99] = [ord(c) for c in struct.pack('<f', self.offsets[0])]  # X0
                    label[100:100+len(self.title)] = [ord(s) for s in self.title]
                    label[244:248] = [ord(c) for c in self.units[0]] + [0]*(4-len(self.units[0]))
                    label[248:252] = [ord(c) for c in self.units[1]] + [0]*(4-len(self.units[1]))
                    label[252:256] = [ord(c) for c in self.units[2]] + [0]*(4-len(self.units[2]))
                # Write label:
                if self.ilabel:
                    f.write(struct.pack('I', 2*256))  # record length, 4 byte format!
                    for element in label:
                        f.write(struct.pack('h', element))  # 2 byte format!
                    f.write(struct.pack('I', 2*256))  # record length!
            # Write picture data:
            for k in range(nlay):
                for j in range(nrow):
                    row = self.data[k, j, :]
                    factor = 8 if self.iform == 3 else 4  # complex numbers need more space!
                    f.write(struct.pack('I', factor*ncol))  # record length, 4 byte format!
                    if self.iform == 0:  # bytes:
                        raise Exception('Byte data is not supported! Use int, float or complex!')
                    elif self.iform == 1:  # int:
                        for element in row:
                            f.write(struct.pack('i', element))  # 4 bytes per data entry!
                    elif self.iform == 2:  # float:
                        for element in row:
                            f.write(struct.pack('f', element))  # 4 bytes per data entry!
                    elif self.iform == 3:  # complex:
                        for element in row:
                            f.write(struct.pack('f', element.real))  # 4 bytes per data entry!
                            f.write(struct.pack('f', element.imag))  # 4 bytes per data entry!
                    f.write(struct.pack('I', factor*ncol))  # record length, 4 byte format!

    @classmethod
    def from_signal(cls, signal):
        # TODO: Docstrings (everywhere)!
        data = signal.data
        assert len(data.shape) <= 3, 'Only up to 3-dimensional datasets can be handled!'
        scales, offsets, units = [1., 1., 1.], [0., 0., 0.], ['', '', '']
        for i in range(len(data.shape)):
            scales[i] = signal.axes_manager[i].scale
            offsets[i] = signal.axes_manager[i].offset
            units[i] = signal.axes_manager[i].units
        for i in range(3 - len(data.shape)):  # Make sure data is 3D!
            data = np.expand_dims(data, axis=0)
        iclass = cls.ICLASS_DICT_INV.get(signal.metadata.Signal.record_by, 6)  # 6: undefined
        if data.dtype.name == 'int8':
            iform = 0  # byte
        elif data.dtype.name in ['int16', 'int32']:
            iform = 1  # int
        elif data.dtype.name in ['float16', 'float32']:
            iform = 2  # float
        elif data.dtype.name == 'complex64':
            iform = 3  #
        else:
            raise TypeError('Data type not understood!')
        arg_dict = {}
        arg_dict['data'] = data
        arg_dict['title'] = signal.metadata.General.as_dictionary().get('title', 'Signal from EMD')
        arg_dict['offsets'] = offsets
        arg_dict['scales'] = scales
        arg_dict['units'] = units
        arg_dict['date'] = strftime('%y-%m-%d %H:%M:%S')
        arg_dict['ICLASS'] = iclass
        arg_dict['IFORM'] = iform
        arg_dict['IVERSN'] = 2  # current standard
        arg_dict['ILABEL'] = 1  # True
        arg_dict['IFORMAT'] = None  # not needed
        arg_dict['IWP'] = 0  # seems standard
        arg_dict['IPLTYP'] = 248  # seems standard
        arg_dict['ICCOLN'] = data.shape[2]//2 + 1
        arg_dict['ICROWN'] = data.shape[1]//2 + 1
        arg_dict['ICLAYN'] = data.shape[0]//2 + 1
        return SemperFormat(arg_dict)

    def to_signal(self):
         # TODO: Docstrings (everywhere)!
        import hyperspy.api as hp
        data = np.squeeze(self.data)  # Reduce unneeded dimensions!
        if self.ICLASS_DICT[self.iclass] == 'spectrum':
            signal = hp.signals.Spectrum(data)
        elif self.ICLASS_DICT[self.iclass] == 'image':
            signal = hp.signals.Image(data)
        else:
            signal = hp.signals.Signal(data)
        for i in range(len(data.shape)):
            signal.axes_manager[i].name = {0: 'x-axis', 1: 'y-axis', 2: 'z-axis'}[i]
            signal.axes_manager[i].scale = self.scales[i]
            signal.axes_manager[i].offset = self.offsets[i]
            signal.axes_manager[i].units = self.units[i]
        signal.metadata.General.title = self.title
        return signal

    def convert_to_abs(self):
        '''Take the absolute of the data. Converts complex to float in the process.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Operates in place!

        '''
        self.data = np.abs(self.data)  # usable for int, float and complex values!
        if self.iform == 3:
            self.iform = 2  # update the format (now it's definitely float)!

    def print_info(self):
        '''Print important flag information of the :class:`.~SemperFormat` object.

        Parameters
        ----------
        None

        Returns
        -------
        None

        '''
        self._log.debug('Calling print_info')
        print '\n------------------------------------------------------'
        print self.title
        print self.date, '\n'
        print 'dimensions: x: {}, y: {}, z: {}'.format(*reversed(self.data.shape))
        print 'scaling:    x: {:.3g}, y: {:.3g}, z: {:.3g}'.format(*self.scales)
        print 'offsets:    x: {:.3g}, y: {:.3g}, z: {:.3g}'.format(*self.offsets)
        print 'units:      x: {}, y: {}, z: {}'.format(*self.units)
        print 'data range:', (self.data.min(), self.data.max()), '\n'
        print 'ICLASS: ', self.ICLASS_DICT[self.iclass]
        print 'IFORM : ', self.IFORM_DICT[self.iform]
        print 'IVERSN: ', self.iversn
        print 'ILABEL: ', self.ilabel == 1
        if self.iformat is not None:
            print 'IFORMAT:', self.iformat
        if self.ilabel:
            print 'IWP   : ', self.iwp
            print 'ICCOLN: ', self.iccoln
            print 'ICROWN: ', self.icrown
            print 'ICLAYN: ', self.iclayn
        print '------------------------------------------------------\n'

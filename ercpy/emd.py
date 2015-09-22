# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 20:58:25 2015

@author: Jan
"""


import re
import h5py

import config

import logging


__all__ = ['EMD']


class EMD(object):

    _log = logging.getLogger(__name__)

    def __init__(self, data={}, user={}, microscope={}, sample={}, comments={}):
        self._log.debug('Calling __init__')
        # Make sure default keys are present in user:
        for key in ['name', 'institution', 'department', 'email']:
            if key not in user:
                user[key] = config.USER_DEFAULTS[key]
        self.user = user
        # Make sure default keys are present in microscope:
        for key in ['name', 'voltage']:
            if key not in microscope:
                microscope[key] = ''
        self.microscope = microscope
        # Make sure default keys are present in sample:
        for key in ['material', 'preparation']:
            if key not in sample:
                sample[key] = ''
        self.sample = sample
        self.comments = comments
        # Make sure the signals are added properly to data:
        self.data = {}
        for name, signal in data.iteritems():
            self.add_signal(name, signal)

    def _write_signal_to_group(self, data_group, signal):
        self._log.debug('Calling _write_signal_to_group')
        # Save data:
        dataset = data_group.create_group(signal.metadata.General.title)
        dataset['data'] = signal.data
        # Iterate over all dimensions:
        for i in range(len(signal.data.shape)):
            # TODO: What about longer dim?
            key = 'dim{}'.format(i+1)
            offset = signal.axes_manager[i].offset
            scale = signal.axes_manager[i].scale
            dim = dataset.create_dataset(key, data=[offset, offset+scale])
            dim.attrs['name'] = signal.axes_manager[i].name
            units = signal.axes_manager[i].units
            dim.attrs['units'] = '[{}]'.format('_'.join(list(units)))
            # TODO: What if units is not one or two characters long? More _? E.g. keV?
        # Write metadata:
        dataset.attrs['emd_group_type'] = 1
        for key, value in signal.metadata.Signal:
            dataset.attrs[key] = value

    def _read_signal_from_group(self, name, group):
        self._log.debug('Calling _read_signal_from_group')
        import hyperspy.hspy as hp
        # Create hyperspy Signal:
        data = group.get('data')[...]
        signal = hp.signals.Signal(data)
        # Set signal properties:
        signal.set_signal_origin = group.attrs.get('signal_origin', '')
        signal.set_signal_type = group.attrs.get('signal_type', '')
        record_by = group.attrs.get('record_by', '')
        # Convert to image or spectrum if necessary:
        if record_by == 'spectrum':
            signal.as_spectrum()
        if record_by == 'image':
            signal.as_image()
        # Iterate over all dimensions:
        for i in range(len(data.shape)):
            dim = group.get('dim{}'.format(i+1))
            signal.axes_manager[i].name = dim.attrs.get('name', '')
            units = re.findall('[^_\W]+', dim.attrs.get('units', ''))
            signal.axes_manager[i].units = ''.join(units)
            try:
                signal.axes_manager[i].scale = dim[1] - dim[0]  # TODO: What about longer dim?
                signal.axes_manager[i].offset = dim[0]
            except Exception:
                print 'Scale could not be calculated!'  # TODO: This happens for strings in dim!
        # Extract metadata:
        metadata = {}
        for key, value in group.attrs.iteritems():
            metadata[key] = value
        # Add signal:
        self.add_signal(name, signal, metadata)

    def add_signal(self, name, signal, metadata={}):
        # TODO: Docstring!
        self._log.debug('Calling add_signal')
        # Save title:
        signal.metadata.General.title = name
        # Save signal metadata:
        signal.metadata.Signal.add_dictionary(metadata)
        # Save global metadata:
        signal.metadata.General.add_node('user')
        signal.metadata.General.user.add_dictionary(self.user)
        signal.metadata.General.add_node('microscope')
        signal.metadata.General.microscope.add_dictionary(self.microscope)
        signal.metadata.General.add_node('sample')
        signal.metadata.General.sample.add_dictionary(self.sample)
        signal.metadata.General.add_node('comments')
        signal.metadata.General.comments.add_dictionary(self.comments)
        # Add signal:
        self.data[name] = signal

    def save_to_emd(self, filename='datacollection.emd'):
        '''Save :class:`~.PhaseMap` data in a file with emd(hdf5)-format.

        Parameters
        ----------
        filename : string, optional
            The name of the emd-file in which to store the phase data.
            The default is '..\output\phasemap.emd'.

        Returns
        -------
        None

        Notes
        -----
        Does not save the unit of the original phase map.

        '''  # TODO: Docstring!
        self._log.debug('Calling save_to_emd')
        # Open file:
        emd_file = h5py.File(filename, 'w')
        # Write version:
        ver_maj, ver_min = config.EMD_VERSION.split('.')
        emd_file.attrs['version_major'] = ver_maj
        emd_file.attrs['version_minor'] = ver_min
        # Write user:
        user_group = emd_file.create_group('user')
        for key, value in self.user.iteritems():
            user_group.attrs[key] = value
        # Write microscope:
        microscope_group = emd_file.create_group('microscope')
        for key, value in self.microscope.iteritems():
            microscope_group.attrs[key] = value
        # Write sample:
        sample_group = emd_file.create_group('sample')
        for key, value in self.sample.iteritems():
            sample_group.attrs[key] = value
        # Write comments:
        comments_group = emd_file.create_group('comments')
        for key, value in self.comments.iteritems():
            comments_group.attrs[key] = value
        # Write data:
        data_group = emd_file.create_group('data')
        for signal in self.data.values():
            self._write_signal_to_group(data_group, signal)
        # Close file and return EMD object:
        emd_file.close()
        return emd

    @classmethod
    def load_from_emd(cls, filename):
        '''Construct :class:`~.PhaseMap` object from NetCDF4-file.

        Parameters
        ----------
        filename : string
            The name of the NetCDF4-file from which to load the data. Standard format is '\*.nc'.

        Returns
        -------
        phase_map: :class:`~.PhaseMap`
            A :class:`~.PhaseMap` object containing the loaded data.

        Notes
        -----
        Does not recover the unit of the original phase map, defaults to `'rad'`.

        '''  # TODO: Docstring!
        cls._log.debug('Calling load_from_ems')
        # Read in file:
        emd_file = h5py.File(filename, 'r')
        # Creat empty EMD instance:
        emd = EMD()
        # Extract user:
        user_group = emd_file.get('user')
        if user_group is not None:
            for key, value in user_group.attrs.iteritems():
                emd.user[key] = value
        # Extract microscope:
        microscope_group = emd_file.get('microscope')
        if microscope_group is not None:
            for key, value in microscope_group.attrs.iteritems():
                emd.microscope[key] = value
        # Extract sample:
        sample_group = emd_file.get('sample')
        if sample_group is not None:
            for key, value in sample_group.attrs.iteritems():
                emd.sample[key] = value
        # Extract comments:
        comments_group = emd_file.get('comments')
        if comments_group is not None:
            for key, value in comments_group.attrs.iteritems():
                emd.comments[key] = value
        # Extract data:
        data_group = emd_file.get('data')
        if data_group is not None:
            for name, group in data_group.iteritems():
                if isinstance(group, h5py.Group):
                    if group.attrs.get('emd_group_type') == 1:
                        emd._read_signal_from_group(name, group)
        # Close file and return EMD object:
        emd_file.close()
        return emd

    def print_info(self):
        self._log.debug('Calling print_info')
        print '\nUser:\n--------------------'
        for key, value in self.user.iteritems():
            print '{}:'.format(key).ljust(15), value
        print '--------------------\n\nMicroscope:\n--------------------'
        for key, value in self.microscope.iteritems():
            print '{}:'.format(key).ljust(15), value
        print '--------------------\n\nSample:\n--------------------'
        for key, value in self.sample.iteritems():
            print '{}:'.format(key).ljust(15), value
        print '--------------------\n\nComments:\n--------------------'
        for key, value in self.comments.iteritems():
            print '{}:'.format(key).ljust(15), value
        print '--------------------\n\nData:\n--------------------'
        for key, value in self.data.iteritems():
            print '{}:'.format(key).ljust(15), value
            print value.metadata
        print '--------------------\n'

# TODO: function to generate subsets of datasets! List as input!

#if __name__ == '__main__':
#    emd = EMD.load_from_emd('C:\Users\Jan\Desktop\EMD\Si3N4_0001_multislice.emd')
#    print emd.data.values()[0].metadata
#    emd.save_to_emd('test.emd')
#    EMD.load_from_emd('test.emd').print_info()

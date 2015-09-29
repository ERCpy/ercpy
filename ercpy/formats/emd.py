# -*- coding: utf-8 -*-
# Copyright 2015 by Forschungszentrum Juelich GmbH
# Author: J. Caron
#
"""This module provides the :class:`~.EMD` class for storing of electron microscopy datasets."""


import re
import h5py

from .. import config

import logging


__all__ = ['EMD']


class EMD(object):

    '''Class for storing electron microscopy datasets.

    The :class:`~.EMD` class can hold an arbitrary amount of datasets in the `data` dictionary.
    These sets are saved as HyperSpy :class:`~hyperspy.signal.Signal` instances. Global metadata
    are saved in four dictionaries (`user`, `microscope`, `sample`, `comments`). To print
    relevant information about the EMD instance use the :func:`~.print_info` function. EMD
    instances can be loaded from and saved to emd-files, an hdf5 standard developed at Lawrence
    Berkeley National Lab (http://emdatasets.lbl.gov/).

    Attributes
    ----------
    data: dictionary
        Dictionary which contains all datasets as :class:`~hyperspy.signal.Signal` instances.
    user: dictionary
        Dictionary which contains user related metadata.
    microscope: dictionary
        Dictionary which contains microscope related metadata.
    sample: dictionary
        Dictionary which contains sample related metadata.
    comments: dictionary
        Dictionary which contains additional commentary metadata.

    '''

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
            name = signal.axes_manager[i].name
            from traits.trait_base import _Undefined
            if type(name) is _Undefined:
                name = ''
            dim.attrs['name'] = name
            units = signal.axes_manager[i].units
            if type(units) is _Undefined:
                units = ''
            else:
                units = '[{}]'.format('_'.join(list(units)))
            dim.attrs['units'] = units
            # TODO: What if units is not one or two characters long? More _? E.g. keV?
        # Write metadata:
        for key, value in signal.metadata.Signal:
            dataset.attrs[key] = value

    def _read_signal_from_group(self, name, group):
        self._log.debug('Calling _read_signal_from_group')
        import hyperspy.api as hp
        # Extract essential data:
        data = group.get('data')[...]
        record_by = group.attrs.get('record_by', '')
        # Create Signal, Image or Spectrum:
        if record_by == 'spectrum':
            signal = hp.signals.Spectrum(data)
        if record_by == 'image':
            signal = hp.signals.Image(data)
        else:
            signal = hp.signals.Signal(data)
        # Set signal properties:
        signal.set_signal_origin = group.attrs.get('signal_origin', '')
        signal.set_signal_type = group.attrs.get('signal_type', '')
        # Iterate over all dimensions:
        for i in range(len(data.shape)):
            dim = group.get('dim{}'.format(i+1))
            signal.axes_manager[i].name = dim.attrs.get('name', '')
            units = re.findall('[^_\W]+', dim.attrs.get('units', ''))
            signal.axes_manager[i].units = ''.join(units)
            try:
                signal.axes_manager[i].scale = dim[1] - dim[0]  # TODO: What about longer dim?
                signal.axes_manager[i].offset = dim[0]
            except (IndexError, TypeError) as e:  # Hyperspy uses defaults (1.0 and 0.0)!
                self._log.warning('Could not calculate scale/offset of axis {}: {}'.format(i, e))
        # Extract metadata:
        metadata = {}
        for key, value in group.attrs.iteritems():
            metadata[key] = value
        # Add signal:
        self.add_signal(name, signal, metadata)

    def add_signal(self, name, signal, metadata={}):
        '''Add a hyperspy signal to the EMD instance and make sure all metadata is present.

        Parameters
        ----------
        name: string
            Name of the (used as a key for the `data` dictionary).
        signal: :class:`~hyperspy.signal.Signal`
            HyperSpy signal which should be added to the EMD instance.
        metadata: dictionary
            Dictionary which holds signal specific metadata which will be added to the signal.

        Returns
        -------
        None

        Notes
        -----
        This is the preferred way to add signals to the EMD instance. Directly adding to the `data`
        dictionary is possible but does not make sure all metadata are correct. This is called
        in the standard constructor on all entries in the `data` dictionary!

        '''
        self._log.debug('Calling add_signal')
        # TODO: Make this ULTRAFLEXIBLE, able to handle unf, dm3, images, txt, all that jazz!
        # TODO: OR have separate add_signal_from_file for that! (which envokes this! jup, better!)
        import hyperspy.api as hp
        if not issubclass(type(signal), hp.signals.Signal):
            try:
                signal = hp.load(signal)  # See if HyperSpy can handle it!
            except Exception as e:
                self._log.error('Could not load signal ({})!'.format(str(e)))
        # Save title:
        signal.metadata.General.title = name
        # Save signal metadata:
        signal.metadata.Signal['emd_group_type'] = 1
        signal.metadata.Signal.add_dictionary(metadata)
        for key in ['name', 'units']:
            if key not in signal.metadata.Signal.as_dictionary():
                signal.metadata.Signal[key] = ''  # make sure at least defaults are present!
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
        '''Save :class:`~.EMD` data in a file with emd(hdf5)-format.

        Parameters
        ----------
        filename : string, optional
            The name of the emd-file in which to store the data.
            The default is 'datacollection.emd'.

        Returns
        -------
        None

        '''
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

    @classmethod
    def load_from_emd(cls, filename):
        '''Construct :class:`~.EMD` object from an emd-file.

        Parameters
        ----------
        filename : string
            The name of the emd-file from which to load the data. Standard format is '*.emd'.

        Returns
        -------
        emd: :class:`~.EMD`
            A :class:`~.EMD` object containing the loaded data.

        '''
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

    # TODO: def export_signal(signal_name, filename):  # Wrappers for several formats!

    def print_info(self):
        '''Print all relevant information about the EMD instance.

        Parameters
        ----------
        None

        Returns
        -------
        None

        '''
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
            print value.metadata.Signal
        print '--------------------\n'

# TODO: function to generate subsets of datasets! List as input!

if __name__ == '__main__':
    emd = EMD.load_from_emd('C:\Users\Jan\Desktop\EMD\Si3N4_0001_multislice.emd')
#    print emd.data.values()[0].metadata
#    emd.save_to_emd('test.emd')
#    EMD.load_from_emd('test.emd').print_info()

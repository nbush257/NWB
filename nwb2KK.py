import pynwb
import re
import os
import glob
import scipy.io.matlab as sio
from datetime import datetime
from dateutil.tz import tzlocal
from pynwb import NWBFile
import numpy as np

#TODO: init_NWB, unpack neuromat,...

def load_neuromat(f):
    ''' quick utility function to load in the matlab data
    with proper organization'''
    dat = sio.loadmat(f,struct_as_record=False,squeeze_me=True)
    dat = dat['data']
    return dat



def init_NWB(p):
    '''
    Initialize an NWB file with appropriate metadata
    :param p: path where all the matlab files exist to be appended to an NWB file

    :return:
    '''
    file_list = glob.glob(os.path.join(p,'*.mat'))
    first_file = os.path.split(file_list[0])[1]
    fname_meta = sanitize_filename(first_file)
    NWBfilename = os.path.join(p,fname_meta['mouse_num'] + '_' + fname_meta['whisker'])
    #TODO init file with relevant metadata

def sanitize_filename(f):
    '''
    takes a filename and makes sure it follows the expected format.
    Then strips all the needed metadata
    :param f: a filename from the matlab conversion script
    :return meta: a dict with metadata stripped from the filename
    '''
    if f[0]!='m':
        raise ValueError('File does not start with "m". Is this a valid mouse file?')

    fileparts = re.split('[_\-]',f)
    if len(fileparts) !=5:
        raise ValueError('Filename does not have the 5 expected parts.')
    if fileparts[-1] != 'tempConvert.mat':
        raise ValueError('filename does not end with tempConvert indicator')
    if len(fileparts[0])!= 7:
        raise ValueError('Mouse number is wrong length')
    if len(fileparts[1])!= 8:
        raise ValueError('Recording date is wrong length')
    if len(fileparts[2])!= 2:
        raise ValueError('Whisker ID is wrong length')
    mouse_num = fileparts[0]
    whisker = fileparts[2]
    trial = int(fileparts[3])

    d_out =  {'whisker':whisker,
            'mouse_num':mouse_num,
            'trial_num':trial}

    return(d_out)


def add_matlab_data(nwb,dat):
    '''
    append data from a matlab file as an experiment in an NWB file

    :param nwb: an NWB fle
    :param dat: a loaded matlab file
    :return: VOID
    '''




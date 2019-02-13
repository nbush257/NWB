import pynwb
from pynwb.ecephys import ElectricalSeries
import pandas as pd
import re
import os
import glob
import scipy.io.matlab as sio
from datetime import datetime,timedelta
from dateutil.tz import tzlocal
from pynwb import NWBFile
import numpy as np
import yaml

#TODO: init_NWB, unpack neuromat,...
#TODO: make abstract loading factories?

def import_NWB_params(yaml_file):
    '''
    Import the parameters of a yaml file that we want to use to init a NWB file
    :param yaml_file: path to the yaml file with the desired parameters
    :return exp_params: a dict of those parameters
    '''
    with open(yaml_file,'r') as fid:
        exp_params = yaml.load(fid)

    return(exp_params)


def load_neuromat(f):
    ''' quick utility function to load in the matlab data
    with proper organization'''
    dat = sio.loadmat(f,struct_as_record=False,squeeze_me=True)
    dat = dat['data']
    return dat


def get_session_starttime(p):
    '''
    Return the start time of all the recordings for a given whisker/day
    :param p: path of the mat files from a session
    :return starttime: a datetime object
    '''
    first_file = get_first_file(p)
    ext = os.splitext(first_file)[-1]
    if ext == '.ddf':
        return(get_ddf_starttime(f))
    else:
        raise ValueError('Unsupported filetype for get_session_starttime: {}'.format(ext))


def get_ddf_starttime(f):

    dat = load_neuromat(os.path.join(p,first_file))
    finfo = dat.fileInfo
    starttime = datetime(finfo.Time_Year,finfo.Time_Month,finfo.Time_Day,finfo.Time_Hour,finfo.Time_Min,finfo.Time_Sec,tzinfo=tzlocal())
    return(starttime)


def get_session_endtime(p):
    '''
    Return the end time of all the recordings for a given whisker/day
    :param p: path of the mat files from a session
    :return endtime: a datetime object
    '''
    last_file = get_last_file(p)
    dat = load_neuromat(os.path.join(p,last_file))
    finfo = dat.fileInfo
    starttime = datetime(finfo.Time_Year,finfo.Time_Month,finfo.Time_Day,finfo.Time_Hour,finfo.Time_Min,finfo.Time_Sec,tzinfo=tzlocal())
    dur = finfo.TimeSpan
    endtime = starttime + timedelta(seconds=dur)
    return(endtime)


def get_first_file(p):
    '''
    Get the first file in a directory. SHould be mat files converted from DDFs
    :param p: directory to search
    :return first_file: the first file in the directory (alphabetical order)
    '''
    file_list = glob.glob(os.path.join(p,'*.mat'))
    first_file = os.path.split(file_list[0])[1]
    return(first_file)


def get_last_file(p):
    '''
    Get the first file in a directory. SHould be mat files converted from DDFs
    :param p: directory to search
    :return first_file: the first file in the directory (alphabetical order)
    '''
    file_list = glob.glob(os.path.join(p,'*.mat'))
    last_file = os.path.split(file_list[-1])[1]
    return(last_file)


def get_subject_info(csv_file,id):
    '''
    read a csv file that contains all the subject data, looking for a specific subject,
    and create an NWB Subject object based on that data.
    :param csv_file: csv file that contains animal data
    :param id: int - id number of the animal you want
    :return: pynwb.file.Subject object populated with the fields from the csv.
    '''
    animal_data = pd.read_csv(csv_file)
    sub_data = animal_data[animal_data['subject_id']==int(id)]
    sub_data.dropna(axis=1,inplace=True)
    sub_data = sub_data.to_dict(orient='records')[0]
    sub_data['subject_id'] = str(sub_data['subject_id'])

    return(pynwb.file.Subject(**sub_data))


def init_NWB(p_load,yaml_file):
    '''
    Initialize an NWB file with appropriate metadata
    :param p: path where all the matlab files exist to be appended to an NWB file
    :return:
    '''
    first_file = get_first_file(p)
    fname_meta = sanitize_filename(first_file)
    NWBfilename = os.path.join(p,fname_meta['mouse_num'] + '_' + fname_meta['whisker'])
    ID = fname_meta['ID']
    starttime = get_session_starttime(p)
    exp_params = import_NWB_params(yaml_file)
    exp_params['subject']=pynwb.file.Subject(subject_id=fname_meta['mouse_num'],genotype='c57bl6')
    nwb_file = NWBFile(exp_params.pop('desc'),ID,
                  session_start_time=starttime,
                  file_create_date=datetime.now(),
                  **exp_params)
    return(nwb_file)


def init_nwb_electrode(electrode_yaml,nwb_file):
    '''
    Creates an electrode object in NWB.
    :param electrode_yaml: a filename of a yaml to use electrode info from
    :return nwb_electrode: an nwb electrode group
    '''
    with open(electrode_yaml,'r') as fid:
        electrode_params = yaml.load(fid)

    device_name = electrode_params['device_name']
    description = electrode_params['description']
    location = electrode_params['location']
    electrode_name = electrode_params['group_name']
    device = nwb_file.create_device(device_name)
    group = nwb_file.create_electrode_group(electrode_name,
                                    location=location,
                                    description=description,
                                    device=device)
    for ii in range(electrode_params['num_sites']):
        nwb_file.add_electrode(x=electrode_params['x'],
                               y=electrode_params['y'],
                               z=electrode_params['z'],
                               imp=electrode_params['imp'],
                               location=location,
                               filtering=electrode_params['filtering'],
                               group=group)


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
    rec_date = fileparts[1]
    whisker = fileparts[2]
    trial = int(fileparts[3])

    d_out =  {'whisker':whisker,
            'mouse_num':mouse_num,
            'trial_num':trial,
            'rec_date':rec_date,
              'ID':'_'.join([mouse_num,whisker,rec_date,'t{:02d}'.format(trial)])}

    return(d_out)


def add_analog_electrode_data(nwb_file,dat):
    '''
    Append the analog electrode data from a matlab file to an nwb file

    :param nwb:
    :param dat:
    :return:
    '''
    electrode_table_region = nwb_file.create_electrode_table_region([0],'this is hardcoded for one electrode') # currently hardcoded for 1 electrode
    ephys_ts = ElectricalSeries('test',
                                dat.analogData.Neural,
                                electrode_table_region,
                                timestamps=dat.time,
                                conversion=1/10000.
                                )
    nwb_file.add_acquisition(ephys_ts)


def append_data(nwb,dat):
    '''
    append data from a matlab file as an experiment in an NWB file

    :param nwb: an NWB fle
    :param dat: a loaded matlab file
    :return: VOID
    '''




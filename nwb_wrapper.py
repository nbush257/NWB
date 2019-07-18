import pynwb
import scipy.stats
from pynwb import NWBHDF5IO
from warnings import warn
from pynwb.ecephys import ElectricalSeries,TimeSeries
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


def overload_ddf(ddf):
    if type(ddf)==str:
        dat = load_neuromat(ddf)
    elif type(ddf) is sio.mio5_params.mat_struct:
        dat = ddf
    return(dat)


def get_ddf_starttime(ddf):
    '''
    Specific start time get for the ddf mat files
    :param ddf: a string to a matlab-ddf file, or a loaded matlab file object
    :return: time of experiment start
    '''
    dat = overload_ddf(ddf)
    finfo = dat.fileInfo
    starttime = datetime(finfo.Time_Year,finfo.Time_Month,finfo.Time_Day,finfo.Time_Hour,finfo.Time_Min,finfo.Time_Sec,finfo.Time_MilliSec,tzinfo=tzlocal())
    starttime = starttime-timedelta(seconds=dat.fileInfo.TimeSpan)

    return(starttime)


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


def init_NWB(f,subject,yaml_file):
    '''
    Initialize an NWB file with appropriate metadata
    :param p: path where all the matlab files exist to be appended to an NWB file
    :return:
    '''
    exp_params = import_NWB_params(yaml_file)

    if type(subject) is dict:
        exp_params['subject']=pynwb.file.Subject(subject_id=subject['ID'],genotype=subject['genotype'])
        ID = subject['ID']
    elif type(subject) is pynwb.file.Subject:
        ID = subject.subject_id
        exp_params['subject'] = subject

    nwb_file = NWBFile(exp_params.pop('desc'),ID,
                  session_start_time=get_ddf_starttime(f),
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


def add_all_analog(ddf,nwb_file):
    '''
    Add all the analog data and information to an NWB as a
    NWB TimeSeries object
    :param ddf: the datawave file or data structure
    :param nwb_file: the NWB file object to write to
    :param starting_time: the start time of this file relative to the first file
    return: None
    '''
    dat = overload_ddf(ddf)
    for label in dat.analogInfo._fieldnames:
        TS = get_analog_TS(dat,label)

        nwb_file.add_acquisition(TS)


def get_analog_TS(dat,label):
    '''
    Converts an individual analog trace from the datawave matlab format to
    the NWB timeseries base class.
    :param dat: a datawave data structure
    :param label: name of the analog data to grab from the ddf data strcut
    :param starting_time: the start time of this file relative to the first file

    :return: TS - a NWB TimeSeries object with the ddf data populated

    '''
    info = dat.analogInfo.__getattribute__(label)
    data = dat.analogData.__getattribute__(label)
    print('Converting {}'.format(label))
    TS = pynwb.base.TimeSeries(label,
                               data=data,
                               unit=info.Units,
                               resolution=info.Resolution,
                               timestamps=dat.time)
    return(TS)


def convert_TS_to_ES(TS,electrode_table_region,gain=10000):
    '''
    This function takes a TimeSeries base class object and converts it to
    and ElectricalSeries subclass object
    :param TS: TimeSeries object
    :param electrode_table_region: the electrode table region that electrode came from
    :param gain: the gain factor of the recording (from the amplifier)
    :return: None
    '''
    ephys_ts = pynwb.ecephys.ElectricalSeries(TS.name,
                                              data=TS.data,
                                              unit=TS.unit,
                                              resolution=TS.Resolution,
                                              timestamps=TS.timestamps,
                                              electrode_table_region=electrode_table_region)
    return(ephys_ts)


def add_all_digital(ddf,nwb_file):
    '''
    Add all the analog data and information to an NWB as a
    NWB TimeSeries object
    :param ddf: the datawave file or data structure
    :param nwb_file: the NWB file object to write to
    :param starting_time: the start time of this file relative to the first file
    return: None
    '''
    dat = overload_ddf(ddf)
    for label in dat.eventData._fieldnames:
        AS = get_digital_AS(dat,label)
        nwb_file.add_acquisition(AS)


def get_digital_AS(dat,label):
    '''
    Converts an individual digital trace from the datawave matlab format to
    the annotation series class
    :param dat: a datawave data structure
    :param label: name of the digital data to grab from the ddf data strcut
    :param starting_time: the start time of this file relative to the first file

    :return: AS - a NWB TimeSeries object with the ddf data populated
    '''
    print('Convertiong {}'.format(label))
    data = dat.eventData.__getattribute__(label)
    if type(data.val) is not np.array:
        data.val = np.array([data.val])
        data.ts = np.array([data.ts])

    AS = pynwb.misc.AnnotationSeries(name=label,
                                     data=data.val,
                                     timestamps=data.ts)
    return(AS)


def get_offset_time(ddf,nwb_file):
    '''
    [deprecated?]
    Utility function to get a given recordings offset from the first recording of the session
    :param ddf: a ddf file or data
    :param nwb_file: the nwb_file to append to
    :return offset_time: float in seconds of the offset between firrst recording start and this recording start
    '''
    # get the time difference between the start of the session and the start of this recording
    dat = overload_ddf(ddf)
    recording_start_time = get_ddf_starttime(dat)
    offset_recording_start_time = recording_start_time - nwb_file.session_start_time
    return(offset_recording_start_time.total_seconds())


def concatenate_recordings(p):
    #[DEPRECATED]
    pass

    '''
    concatenate all the recordings in a given path
    :param p: Path where all the ddf mat files live
    :return cat_dict : a dictionary of the concatenated data
    '''
    recording_start_times = []
    recording_lengths = []
    recording_start_indices = [0]
    nframes = []
    trial_nums = []
    for ii,f in enumerate(glob.glob(os.path.join(p,'*.mat'))):
        # ignore calibration files
        if re.search('calib',f) is not None:
            print('Skipping {}'.format(os.path.split(f)[-1]))
            continue
        print('Loading {}'.format(os.path.split(f)[-1]))
        fname_meta = sanitize_AWAKE_filename(f)
        trial_nums.append(fname_meta['trial_num'])
        dat = overload_ddf(f)
        if ii==0:
            # if this is the first file in the session, initialize the time basis
            ts = dat.time
            recording_start_times.append(ts[0])
            starttime = get_ddf_starttime(dat)
            recording_lengths.append(dat.fileInfo.TimeSpan)

            # create signals
            analog={}
            for field in dat.analogData._fieldnames:
                analog[field] = dat.analogData.__getattribute__(field)

        else:
            for field in dat.analogData._fieldnames:
                analog[field] = np.concatenate([analog[field],dat.analogData.__getattribute__(field)])

            # get the start time of this recording, calculate the offset from the first recording, and append to the end
            exp_time = get_ddf_starttime(dat)
            offset = exp_time-starttime
            offset_time = dat.time+offset.total_seconds()
            recording_start_times.append(offset.total_seconds())
            recording_start_indices.append(len(ts))
            recording_lengths.append(dat.fileInfo.TimeSpan)
            # get the camera trigger for this recording and offset according to first recording
            #exp_idx = trigger_to_idx(dat.analogData.Cam_trig)
            #exp_offset_idx = exp_idx+len(ts)
            #exp_offset_frametimes = dat.time[exp_idx] + offset.total_seconds()

            # if the frametimes are not very regular, then the trigger probably failed
            #if np.var(np.diff(exp_offset_frametimes[1:]))>1e-5:
                #exp_offset_idx = np.array([],dtype='int')
                #exp_offset_frametimes = np.array([])
                #warn('Variance of frametimes is high. Probably not the actual camera trigger.\n Deleting all frames in {}'.format(f))
            #nframes.append(len(exp_offset_idx))

            # append offset neural and triggers to the full dataset
            ts = np.concatenate([ts,offset_time])
            #frame_idx = np.concatenate([frame_idx,exp_offset_idx])
            #frametimes = np.concatenate([frametimes,exp_offset_frametimes])

    cat_dict = {'ts': ts,
                'analogData': analog,
                'recording_start_times':np.array(recording_start_times),
                'recording_start_indices':np.array(recording_start_indices),
                'whisker':fname_meta['whisker'],
                'trial_nums':trial_nums,
                'subject_id':fname_meta['mouse_num'],
                'rec_date':fname_meta['rec_date'],
                'recording_lengths':recording_lengths
                }
    return(cat_dict)


def convert_NWB(f,exp_yaml,subject=None,NWBfilename=None):
    '''
    convert a Matlab DDF to a NWB formatted file
    :param f: DDF file to convert
    :param exp_yaml: a yaml which holds metadata for the experiment
    :param subject: a pynwb subject object
    :param NWBfilename: the output filename to convert to. If not provided, simply changes the extension
    :return:
    '''
    p = os.path.split(f)[0]

    # If no output name is supplied. simply replace the extension with nwb
    if NWBfilename is None:
        NWBfilename = os.path.splitext(f)[0]+'.nwb'

    # If no subject info is given, ask user for ID and genotype
    if subject is None:
        subject={}
        subject['ID'] = input('Enter the subject ID number')
        subject['genotype'] = input('Enter the subject genotype (default=c57bl6)')
        if len(subject['genotype'])==0:
            subject['genotype'] = 'c57bl6'


    # Load matlab DDF and add all entries
    dat = overload_ddf(f)
    nwb_file = init_NWB(f,subject,exp_yaml)
    add_all_analog(dat,nwb_file)
    add_all_digital(dat,nwb_file)

    # write the file
    print('writing NWB_file to\n\t{}...'.format(NWBfilename))
    with NWBHDF5IO('{}'.format(NWBfilename), 'w') as io:
        io.write(nwb_file)
    print('Done!')



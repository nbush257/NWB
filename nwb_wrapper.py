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


def get_session_starttime(p):
    '''
    Return the start time of all the recordings for a given whisker/day
    :param p: path of the mat files from a session
    :return starttime: a datetime object
    '''
    first_file = get_first_file(p)
    ext = os.path.splitext(first_file)[-1]
    if ext == '.mat':
        return(get_ddf_starttime(first_file))
    else:
        raise ValueError('Unsupported filetype for get_session_starttime: {}'.format(ext))


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


def get_first_file(p,fullfile_tgl=True):
    '''
    Get the first file in a directory. SHould be mat files converted from DDFs
    :param p: directory to search
    :return first_file: the first file in the directory (alphabetical order)
    '''
    file_list = glob.glob(os.path.join(p,'*.mat'))
    if not fullfile_tgl:
        first_file = os.path.split(file_list[0])[1]
    else:
        first_file = file_list[0]
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


def init_NWB(p,yaml_file):
    '''
    Initialize an NWB file with appropriate metadata
    :param p: path where all the matlab files exist to be appended to an NWB file
    :return:
    '''
    first_file = get_first_file(p)
    fname_meta = sanitize_AWAKE_filename(first_file)
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


def sanitize_AWAKE_filename(f):
    '''
    takes a filename and makes sure it follows the expected format.
    Then strips all the needed metadata
    :param f: a filename from the matlab conversion script
    :return meta: a dict with metadata stripped from the filename
    '''
    f = os.path.split(f)[-1]
    if f[0]!='m':
        raise ValueError('File does not start with "m". Is this a valid mouse file?')

    fileparts = re.split('[_\-]',f)
    if len(fileparts) !=5:
        raise ValueError('Filename does not have the 5 expected parts.')
    if fileparts[-1] != 'tempConvert.mat':
        raise ValueError('filename does not end with tempConvert indicator')
    if len(fileparts[0])!= 7:
        raise ValueError('Mouse number is wrong length')
    if len(fileparts[2])!= 8:
        raise ValueError('Recording date is {}, and is of wrong length'.format(fileparts[1]))
    if len(fileparts[1])!= 2:
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


def add_all_analog(ddf,nwb_file,starting_time):
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
        TS = get_analog_TS(dat,label,starting_time)

    nwb_file.add_acquisition(TS)


def get_analog_TS(dat,label,starting_time):
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
                               timestamps=dat.time,
                               starting_time=starting_time)
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
                                              starting_time=TS.starting_time,
                                              electrode_table_region=electrode_table_region)
    return(ephys_ts)


def add_all_digital(ddf,nwb_file,starting_time):
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
        AS = get_digital_AS(dat,label,starting_time)

    nwb_file.add_acquisition(AS)


def get_digital_AS(dat,label,starting_time):
    '''
    Converts an individual digital trace from the datawave matlab format to
    the annotation series class
    :param dat: a datawave data structure
    :param label: name of the digital data to grab from the ddf data strcut
    :param starting_time: the start time of this file relative to the first file

    :return: AS - a NWB TimeSeries object with the ddf data populated
    '''
    data = dat.eventData.__getattribute__(label)

    AS = pynwb.misc.AnnotationSeries(name=label,
                                     data=data.val,
                                     timestamps=data.ts+starting_time,
                                     )
    return(AS)


def mark_positions(motor,nwb_file):
    '''
    Add the times during which the motor was in a given position
    i.e., exclude times when motor is moving.
    :param motor: The digital read from the motor
    :param nwb_file: the NWB file object to write to
    :return:
    '''
    temp = np.concatenate([[0.],motor.timestamps])
    position_times = np.reshape(temp[:-1],[-1,2]) # marks the motor motion time
    for ii,(onset,offset) in enumerate(position_times):
        nwb_file.add_epoch(onset,offset,'pos_{0:03d}'.format(ii))


def add_neural_analog(f,nwb_file):
    if os.path.splitext(f)[-1]=='.mat':
        get_AWAKE_neural(f,nwb_file)
    else:
        raise ValueError('Filetype not supportd to add analog signal: {}'.format(os.path.splitext(f)[-1]))


def get_AWAKE_neural(ddf,nwb_file,gain=10000):
    '''
    extract the neural data from a matlab ddf. Assumes one channel of recording,
    and that the one channel has a specific format in the matlab ddf.
    :param f: a matlab ddf
    :param conversion: the gain of the amplifier
    :return ephys_ts: a NWB ElectricalSeries object
    '''
    dat = overload_ddf(ddf)
    electrode_table_region = nwb_file.create_electrode_table_region([0],
                                'this is hardcoded for one electrode') # currently hardcoded for 1 electrode

    # get the time difference between the start of the session and the start of this recording

    offset_time = get_offset_time(dat,nwb_file)
    # create the ephys object and add it to the nwb file
    ephys_ts = ElectricalSeries('channel_0',
                                dat.analogData.Neural,
                                electrode_table_region,
                                timestamps=dat.time+offset_time,
                                conversion=1/float(gain)
                                )
    nwb_file.add_acquisition(ephys_ts)
    return(ephys_ts)


def add_camera_trigger(ddf,nwb_file):
    '''
    use this to transform the camera trigger analog input into events
    :param ddf:
    :param nwb_file:
    :return:
    '''

    dat = overload_ddf(ddf)
    offset_time = get_offset_time(dat,nwb_file)
    idx = trigger_to_idx(dat.analogData.Cam_trig)
    ts = dat.time[idx]
    frametimes = pynwb.ecephys.TimeSeries('Camera Frame onsets',
                                          idx,
                                          'indices',
                                          timestamps=ts+offset_time)
    nwb_file.add_acquisition(frametimes)
    return(frametimes)


def trigger_to_idx(trigger,thresh=100.,sign='+'):
    '''
    Takes an analog signal of a camera trigger and finds the
    when the trigger crosses.
    :param trigger:
    :param thresh:
    :return idx: indices in the analog signal in which the threshold was croseed
    '''
    # get the initial value of the trigger, along with the 10th percentile lagrest and smallest values
    initial = np.median(trigger[:100])
    biggest = np.quantile(trigger,.9)
    smallest = np.quantile(trigger,.1)

    if (biggest-smallest) <1000:
        raise ValueError('Difference between big and small values in camera trigger is less than 1000')
    thresh = np.mean([biggest,smallest])

    above = trigger>thresh
    d = np.diff(above.astype('double'))
    first = np.where(np.abs(d)>.5)[0][0]
    if d[first]>0:
        idx =np.where(d>.5)[0]
    else:
        idx =np.where(d<-.5)[0]

    frame_interval = np.mean(np.diff(idx))/40000. #calculate the frame interval in seconds assuming a 40kHz sample
    print('='*15)
    print('Frame rate of camera is calculated to be {:.2f} fps\n'.format(frame_interval**-1))

    return(idx)



def get_offset_time(ddf,nwb_file):
    '''
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


def concatenate_AWAKE_recordings(p):

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
            ndata = dat.analogData.Neural
            ts = dat.time
            recording_start_times.append(ts[0])
            starttime = get_ddf_starttime(dat)
            frame_idx = trigger_to_idx(dat.analogData.Cam_trig)
            frametimes = dat.time[frame_idx]
            recording_lengths.append(dat.fileInfo.TimeSpan)
            # if the frametimes are not very regular, then the trigger probably failed
            if np.var(np.diff(frametimes[1:]))>1e-5:
                frame_idx = np.array([],dtype='int')
                frametimes = np.array([])
                warn('Variance of frametimes is high. Probably not the actual camera trigger.\n Deleting all frames in {}'.format(f))
            nframes.append(len(frame_idx))
        else:
            ndata = np.concatenate([ndata,dat.analogData.Neural])

            # get the start time of this recording, calculate the offset from the first recording, and append to the end
            exp_time = get_ddf_starttime(dat)
            offset = exp_time-starttime
            offset_time = dat.time+offset.total_seconds()
            recording_start_times.append(offset.total_seconds())
            recording_start_indices.append(len(ts))
            recording_lengths.append(dat.fileInfo.TimeSpan)
            # get the camera trigger for this recording and offset according to first recording
            exp_idx = trigger_to_idx(dat.analogData.Cam_trig)
            exp_offset_idx = exp_idx+len(ts)
            exp_offset_frametimes = dat.time[exp_idx] + offset.total_seconds()

            # if the frametimes are not very regular, then the trigger probably failed
            if np.var(np.diff(exp_offset_frametimes[1:]))>1e-5:
                exp_offset_idx = np.array([],dtype='int')
                exp_offset_frametimes = np.array([])
                warn('Variance of frametimes is high. Probably not the actual camera trigger.\n Deleting all frames in {}'.format(f))
            nframes.append(len(exp_offset_idx))

            # append offset neural and triggers to the full dataset
            ts = np.concatenate([ts,offset_time])
            frame_idx = np.concatenate([frame_idx,exp_offset_idx])
            frametimes = np.concatenate([frametimes,exp_offset_frametimes])

    cat_dict = {'ts': ts,
                'frame_idx': frame_idx,
                'frame_times': frametimes,
                'neural': ndata,
                'recording_start_times':np.array(recording_start_times),
                'recording_start_indices':np.array(recording_start_indices),
                'nframes':nframes,
                'whisker':fname_meta['whisker'],
                'trial_nums':trial_nums,
                'subject_id':fname_meta['mouse_num'],
                'rec_date':fname_meta['rec_date'],
                'recording_lengths':recording_lengths
                }
    return(cat_dict)


def add_trial_data(nwb_file,cat_dict):
    '''
    Add information on when trials occured into the nwb file
    :param nwb_file: file to add trial information to
    :param cat_dict: dict returned when you concatenated the ddfs
    :return: None
    '''
    #TODO: import trial information from external source
    for trial_start,dur in zip(cat_dict['recording_start_times'],cat_dict['recording_lengths']):
        nwb_file.add_trial(start_time=trial_start,stop_time=trial_start+dur)



def write_AWAKE_NWB(p,exp_yaml,electrode_yaml,gain=10000.):
    '''
    Create an NWB file for the all the recordings of a particular session for the AWAKE data
    :param p: path where all the ddf-mat files live
    :param exp_yaml: a yaml which holds metadata for the experiment
    :return:
    '''
    first_file = get_first_file(p)
    fname_meta = sanitize_AWAKE_filename(first_file)
    NWBfilename = os.path.join(p,fname_meta['mouse_num'] + '_' + fname_meta['rec_date'] + '_' + fname_meta['whisker']+'.nwb')
    cat_dict = concatenate_AWAKE_recordings(p)
    nwb_file = init_NWB(p,exp_yaml)
    init_nwb_electrode(electrode_yaml,nwb_file)
    electrode_table_region = nwb_file.create_electrode_table_region([0],
                                                                    'this is hardcoded for one electrode') # currently hardcoded for 1 electrode
    ephys_ts = ElectricalSeries('channel_0',
                                cat_dict['neural'],
                                electrode_table_region,
                                timestamps=cat_dict['ts'],
                                conversion=1./float(gain))
    nwb_file.add_acquisition(ephys_ts)
    cam_trig = TimeSeries('Camera frame onsets',
                          cat_dict['frame_idx'],
                          'indices',
                          timestamps=cat_dict['frame_times'])
    nwb_file.add_acquisition(cam_trig)
    add_trial_data(nwb_file,cat_dict)
    # write the file
    print('writing NWB_file to\n\t{}...'.format(NWBfilename))
    with NWBHDF5IO('{}'.format(NWBfilename), 'w') as io:
        io.write(nwb_file)
    print('Done!')



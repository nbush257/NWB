import numpy as np
from tqmd import tqdm
from nwb_wrapper import *
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
    if len(fileparts[1])!= 8:
        raise ValueError('Recording date is {}, and is of wrong length'.format(fileparts[1]))
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


def mark_positions(motor,nwb_file):
    '''
    Add the times during which the motor was in a given position
    i.e., exclude times when motor is moving.
    :param motor: The digital read from the motor (A pynwb acquisition
    :param nwb_file: the NWB file object to write to
    :return:
    '''
    mdata = motor.data.astype('float').astype('int')
    idx = np.logical_or(mdata == 0, mdata == 1)
    ts = motor.timestamps.value[idx]

    temp = np.concatenate([[0.],ts])
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

def ts_to_index(spike_times,sample_times):
    spike_frame = np.empty_like(spike_times)
    for ii,spike in tqdm(enumerate(spike_times)):
        spike_frame[ii] = np.argmax(spike<sample_times)
    return(spike_frame)


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

def create_AWAKW():

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

import pynwb
import subprocess
import matplotlib.pyplot as plt
import sys
import numpy as np
import os
from pynwb import NWBHDF5IO
import h5py
from optparse import OptionParser
from scipy.signal import butter


def write_raw_for_kk(nwb_filename,save_dir=None,detect_direction=None):
    '''
    Write a raw file the Klusta kwik can work on
    :return:
    '''
    if save_dir is None:
        save_dir = os.path.split(nwb_filename)[0]
    if not os.path.isdir(save_dir):
        print('='*10)
        print('Save directory was not found, making directory at:\n\t{}'.format(save_dir))

        os.makedirs(save_dir)

    basename = os.path.splitext(os.path.split(nwb_filename)[-1])[0]
    out_filename = os.path.join(save_dir,basename+'.dat')
    if os.path.isfile(out_filename):
        print('\n'+'='*15)
        print('File {} already exists\n'.format(out_filename))
        overwrite_flag = input('Overwrite(Y/n)')
        if overwrite_flag.lower() != 'y':
            print('Exiting without writing...')
            return(0)


    with NWBHDF5IO(nwb_filename,'r') as io:
        nwbfile = io.read()
        ephys_ts = nwbfile.get_acquisition('single channel of neural data across recordings')
        if detect_direction is None:
            detect_direction = manual_detect_direction(ephys_ts)
        raw_array = np.array(ephys_ts.data,dtype = 'float32')
        print('\n'+'='*15)
        print('Writing raw data to {}'.format(out_filename))
        print('\n'+'='*15)
        raw_array.tofile(out_filename)
        print('done')

    return(detect_direction)


def manual_detect_direction(ephys_ts):
    '''
    Get user input to get spike direction
    :param ephys_ts:
    :return: string indicating detect direction
    '''
    #TODO: Plotting is pretty restrictive. Could make it better if need be
    plt.plot(ephys_ts.data[:1000000],'k')
    plt.title('Click above or below zero to indicate spike direction')
    detect_direction = plt.ginput(1)[0][1]
    plt.close('all')
    if detect_direction>0.:
        detect_direction='positive'
    elif detect_direction<0.:
        detect_direction='negative'
    return(detect_direction)

def kk2nwb_units(kwik_filename,nwb_filename):
    '''
    Use this function to read sorted spikes from a kwik file and append them to the nwb file
    :return None:
    '''

    # open the nwb file for reading and appending
    io = NWBHDF5IO(nwb_filename,mode='a')
    nwbfile = io.read()
    time_vec = nwbfile.acquisition['single channel of neural data across recordings'].timestamps.value
    neural = nwbfile.acquisition['single channel of neural data across recordings'].data

    # access the kwik sort data
    kwik = h5py.File(kwik_filename)
    spikesamps = kwik['channel_groups/0/spikes/time_samples'].value
    spike_cluster_id = kwik['channel_groups/0/spikes/clusters/main'].value
    clusters = kwik['channel_groups/0/clusters/main']

    # extract the cluster ids that were labelled good
    good = []
    num_spikes=[]
    for cluster in clusters:
        if clusters[cluster].attrs['cluster_group']==2:
            good.append(int(cluster))
            num_spikes.append(len(np.where(spike_cluster_id==int(cluster))[0]))

    # sort good descending by number of spikes in the cluster
    good = [x for _,x in sorted(zip(num_spikes,good))][::-1]
    num_spikes = sorted(num_spikes)[::-1]

    print('\nFound {} good clusters\n'.format(len(good)))
    for clu,num in zip(good,num_spikes):
        print('\tcluster {}:\tnum_spikes: {}\n'.format(clu,num))
    # get spike_id for spikes labelled in each good cluster and add them to units in the NWB file
    for ii,clu in enumerate(good):
        spike_idx = spikesamps[spike_cluster_id==clu]
        wvfm_m,wvfm_sd = get_waveforms(spike_idx,neural)[1:]

        nwbfile.add_unit(id=ii,electrodes=[0],
                         spike_times=time_vec[spikesamps],
                         waveform_mean=wvfm_m,
                         waveform_sd=wvfm_sd)

    print('Writing unit data to {}...'.format(nwb_filename))
    # write and close the nwb file
    io.write(nwbfile)
    io.close()
    print('done!')
    return(0)

#TODO: remove klustakwik leftovers


def get_neural(raw_dat_file,dat_type='float32'):
    '''
    Extract the raw neural data from the linked dat file
    :param raw_dat_file:
    :return:
    '''
    raw = np.fromfile(raw_dat_file,dat_type)
    return(raw)


def get_waveforms(spike_idx,neural,pre=15,post=15):
    '''
    Extract the mean and standard deviation of spikes that occur at given indices in the
    given neural data
    :param spike_idx: indices of the spikes detection
    :param neural: raw or filtered neural data from which to extract the waveshapes
    :param pre: window before detection to extract
    :param post: window after detection to extract
    :return: waveforms,waveform_mean,waveform_sd
    '''
    pre = np.asarray(pre,dtype='uint64')
    post = np.asarray(post,dtype='uint64')
    waveforms = np.empty([pre+post,len(spike_idx)])
    for ii,idx in enumerate(spike_idx):
        waveforms[:,ii] = neural[idx-pre:idx+post]
    return(waveforms,np.mean(waveforms,axis=1),np.std(waveforms,axis=1))


def FHC_SE_prb_file(prb_filename):
    """
    write a probe file that has the parameters for a single electrode recording
    :param prb_filename: filename for the probe file
    :return:
    """
    fid = open(prb_filename,'w')
    fid.write("channel_groups={0: {'channels': [0], 'graph': [(0, 0)], 'geometry': {0: (0, 0)}}}")
    fid.close()
    return(0)

def save_param_file(nwb_filename,detect='positive'):
    '''
    saves a parameter file to the location of the nwb file to sort.
    currently hardcoded some default parameters, while reading in some from the file
    :return:
    '''
    if detect not in ['positive','negative']:
        raise ValueError('Detect must be "positive" or "negative"')

    # TODO: make SR calculation better
    io = NWBHDF5IO(nwb_filename,'r')
    nwbfile = io.read()
    param_filename = os.path.splitext(nwb_filename)[0]+'.prm'
    prb_filename = os.path.join(os.path.split(nwb_filename)[0],'se.prb')
    print('Probe filename is {}'.format(prb_filename))
    if not os.path.isfile(prb_filename):
        FHC_SE_prb_file(prb_filename)

    print('Sample Rate conversion is kinda hacky\n')

    time_vec = nwbfile.acquisition['single channel of neural data across recordings'].timestamps.value
    sr = np.mean(np.diff(time_vec[:1000])) ** -1
    gain = nwbfile.acquisition['single channel of neural data across recordings'].conversion

    print('Writing parameter file to: {}'.format(param_filename))
    fid = open(param_filename,'w')
    experiment_name = os.path.splitext(os.path.split(nwb_filename)[1])[0]
    fid.write('experiment_name = "{}"\n'.format(experiment_name))
    fid.write('prb_file="{}"\n'.format(os.path.split(prb_filename)[1]))
    fid.write('traces=dict(raw_data_files=["{}"],voltage_gain={},sample_rate={},n_channels=1,dtype="{}")\n'.format(
        experiment_name+'.dat',
        gain,
        sr,
        'float32'))

    #spikedetekt params
    fid.write('spikedetekt={}\n')
    fid.write("spikedetekt['filter_low']={}\n".format(300.))
    fid.write("spikedetekt['filter_high_factor']={}\n".format(0.95*.5))
    fid.write("spikedetekt['butter_order']={}\n".format(3))
    fid.write("spikedetekt['chunk_size_seconds']={}\n".format(1))
    fid.write("spikedetekt['chunk_overlap_seconds']={}\n".format(0.015))
    fid.write("spikedetekt['n_excerpts']={}\n".format(50))
    fid.write("spikedetekt['excerpt_size_seconds']={}\n".format(1))
    fid.write("spikedetekt['threshold_strong_std_factor']={}\n".format(6.5))
    fid.write("spikedetekt['threshold_weak_std_factor']={}\n".format(4.))
    fid.write("spikedetekt['detect_spikes']='{}'\n".format(detect))
    fid.write("spikedetekt['connected_component_join_size']={}\n".format(1))
    fid.write("spikedetekt['extract_s_before']={}\n".format(15))
    fid.write("spikedetekt['extract_s_after']={}\n".format(15))
    fid.write("spikedetekt['n_features_per_channel']={}\n".format(5))
    fid.write("spikedetekt['pca_n_waveforms_max']={}\n".format(10000))

    # Klustakwik Params
    fid.write("klustakwik2 = {}\n")
    fid.write("klustakwik2['num_starting_clusters'] = {}\n".format(25))

    fid.close()
if __name__=='__main__':
    usage = "Usage: %prog [options] nwb_filename"
    parser = OptionParser(usage=usage)
    parser.add_option('-s','--savedir',
                      type='string',
                      dest='save_dir',
                      action='store',
                      default=None,
                      help='Directory where we save the raw neural data extracted from the nwb file')
    parser.add_option('-k','--kwik_filename',
                      type='string',
                      dest='kwik_filename',
                      action='store',
                      default=None,
                      help='Filename of the sorted kwik file. If this argument is passed, then appends to an nwb file')

    parser.add_option('-d','--detect_direction',
                      type='string',
                      dest='detect',
                      action='store',
                      default=None,
                      help='[deprecated] Direction of the spike to detect')

    options,args = parser.parse_args()
    if (options.save_dir is not None) and (options.kwik_filename is not None):
        raise ValueError('Both a kwik file and a save directory was passed in command line.\n Both options are not used simultaneously')
    # print('Options:\n',options)
    # print('NWB_file:',args)
    if len(os.path.split(args[0])[0])==0:
        p = os.getcwd()
        args[0] = os.path.join(p,args[0])

    print(args)
    if options.kwik_filename is not None:
        kk2nwb_units(options.kwik_filename,args[0])
    else:
        detect_direction = write_raw_for_kk(args[0],save_dir=options.save_dir,detect_direction=options.detect)
        save_param_file(args[0],detect=detect_direction)
        prm_file = os.path.splitext(args[0])[0]
        print('Running klusta on {}.prm'.format(prm_file))
        os.system('deactivate & activate klusta & klusta {}.prm'.format(prm_file))











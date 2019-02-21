import pynwb
import sys
import numpy as np
import os
from pynwb import NWBHDF5IO
import h5py
from optparse import OptionParser
from scipy.signal import butter


def write_raw_for_kk(nwb_filename,save_dir=None):
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
        raw_array = np.array(ephys_ts.data,dtype = 'float32')
        print('\n'+'='*15)
        print('Writing raw data to {}'.format(out_filename))
        print('\n'+'='*15)
        raw_array.tofile(out_filename)
        print('done')

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


# TODO: Run klusta kwik from python call?
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

    options,args = parser.parse_args()
    if (options.save_dir is not None) and (options.kwik_filename is not None):
        raise ValueError('Both a kwik file and a save directory was passed in command line.\n Both options are not used simultaneously')
    # print('Options:\n',options)
    # print('NWB_file:',args)

    if options.kwik_filename is not None:
        kk2nwb_units(options.kwik_filename,args[0])
    else:
        write_raw_for_kk(args[0],save_dir=options.save_dir)










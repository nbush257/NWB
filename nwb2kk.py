import pynwb
import os
from pynwb import NWBHDF5IO
import h5py


def write_raw_for_kk(f):
    '''
    Write a raw file the Klusta kwik can work on
    :return:
    '''
    nwbfile = NWBHDF5IO(f,'r').read()
    ephys_ts = nwbfile.get_acquisition('single channel of neural data across recordings')
    raw_array = np.array(ephys_ts.data,dtype = 'float32')
    raw_array.tofile(r'C:\Users\guru\Desktop\test\test.dat')



def kk2nwb_units(p,prefix,nwbfile):
    '''
    Use this function to read sorted spikes from a kwik file.
    :return:
    '''
    kwik = h5py.File(os.path.join(p,prefix+'.kwik'))
    kwx = h5py.File(os.path.join(p,prefix+'.kwx'))
    spikesamps = kwik['channel_groups/0/spikes/time_samples'].value
    spike_cluster_id = kwik['channel_groups/0/spikes/clusters/main'].value
    clusters = kwik['channel_groups/0/clusters/main']
    good = [] # get cluster_id for spikes labelled good
    for cluster in clusters:
        if clusters[cluster].attrs['cluster_group']==2:
            good.append(int(cluster))

    # get spike_id for spikes labelled in each good cluster
    for ii in good:
        idx = spike_cluster_id==ii









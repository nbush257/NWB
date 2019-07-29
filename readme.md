## This is a readme for converting ddf files to nwb format and the sorting those files.
### Be sure to activate the NWB conda environment
From the command line: `activate nwb`
### Step by step workflow:
1. Place all the ddf files you want to work with into a directory.
1. cd into that directory in Matlab and run `ddf2mat.m`.  
This will create convert the DDF proprietary data format into a Matlab struct. This file is meant to be a temporary conversion format. 
You should now have a new directory with a list of mat files. We will call this "data_path" 
1. activate the NWB environment with `activate nwb`
1. Run `python batch_mat2nwb.py "data_path" -e "exp_yaml" -c -i "`  
Where `exp_yaml` is a text file that contains the experimental metadata  
This converts the matlab ddfs into individual nwb files, concatenates them, and ignores calibration files. You will be asked at the end if you want to remove the old individual files.  
^*This is suggested*
1. Run the `toAwake.ipynb` notebook being sure to edit the needed filenames and parameters
1. Run `python nwb2kk.py [concatenated_nwb.nwb]`   
This sorts the data
1. Once the data is sorted, run `python nwb2kk.py [concatenated_nwb.nwb] -k [sorted_kwik_file.kwik]`  
This adds the sorted unit data to the nwb file.

#### Now you have the sorted unit in the concatenated NWB file.
In python, to access the data, use the following functions:
```python
import pynwb
io = pynwb.NWBHDF5IO(nwb_file,'r')
data = io.read()
neural_data = data.get_acquisition('Neural').data
neural_ts = data.get_acquisition('Neural').timestamps

unit_num = 0
num_units = len(data.units)
spikesamps = data.units['spike_times'][unit_num]
avg_waveform = data.units['waveform_mean'][unit_num]

frame_samps = data.get_acquisition('Trigger').data.value
frame_time = data.get_acquisition('Trigger').timestamps.value


```






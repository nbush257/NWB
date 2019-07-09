This is a readme for converting ddf files to nwb format and the sorting those files.

First, place all the ddf files you want to work with into a directory.
Second, cd into that directory in Matlab and run "ddf2mat.m". This will create convert the DDF proprietary data format into a Matlab struct. This file is meant to be a temporary conversion format. 
You should now have a new directory with a list of mat files. We will call this "data_path" 
We next want to convert those many .mat files into one NWB file. This is done by running mat2nwb which calls functions from nwb_wrapper.






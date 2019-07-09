from nwb_wrapper import *
from optparse import OptionParser

if __name__=='__main__':
    usage = "Usage: %prog [options] path exp_yaml electrode_yaml"
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
    print(args)
#
#TODO: use concatenated awake recordings to create an NWB file
#TODO: make the NWB file readable by KlustaKwik


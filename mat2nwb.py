from nwb_wrapper import *
from optparse import OptionParser

if __name__=='__main__':
    usage = "Usage: %prog path -e exp_yaml -l electrode_yaml"
    parser = OptionParser(usage=usage)
    parser.add_option('-e','--experiment',
                      type='string',
                      dest='exp_yaml',
                      action='store',
                      default=None,
                      help='Yaml file with the experiment paramters')
    parser.add_option('-l','--electrode',
                      type='string',
                      dest='electrode_yaml',
                      action='store',
                      default=None,
                      help='Yaml file with the elcetrode paramters')


    options,args = parser.parse_args()

    print(options)
    print(args)
    write_AWAKE_NWB(args[0],options.exp_yaml,options.electrode_yaml)
#
#TODO: use concatenated awake recordings to create an NWB file
#TODO: make the NWB file readable by KlustaKwik


from nwb_wrapper import *
from optparse import OptionParser

if __name__=='__main__':
    usage = "Usage: %prog path -e exp_yaml -s subject -i ID "
    parser = OptionParser(usage=usage)
    parser.add_option('-e','--experiment',
                      type='string',
                      dest='exp_yaml',
                      action='store',
                      default=None,
                      help='Yaml file with the experiment paramters')
    parser.add_option('-s','--subject',
                      type='string',
                      dest='subject',
                      action='store',
                      default=None,
                      help='Subject ID')


    options,args = parser.parse_args()

    print(options)
    print(args)

    if options.subject is None:
    # If no subject info is given, ask user for ID and genotype
        subject = input('Enter the subject ID number')

    file_list = glob.glob(os.path.join(args[0],'*.mat'))
    print(*file_list)
    print(subject)
    for f in file_list:
        convert_NWB(f=f,
                    exp_yaml=options.exp_yaml,
                    ID=subject)
#
#TODO: use concatenated awake recordings to create an NWB file
#TODO: make the NWB file readable by KlustaKwik


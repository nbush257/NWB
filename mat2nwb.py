from nwb_wrapper import *
from optparse import OptionParser

if __name__=='__main__':
    usage = "Usage: %prog matfile -e exp_yaml -s subject -i ID -o outfilename"
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
                      help='A csv with the subject information')
    parser.add_option('-i','--id',
                      type='string',
                      dest='subject_id',
                      action='store',
                      default=None,
                      help='The ID number of the subject to grab from the csv')
    parser.add_option('-o','--out',
                      type='string',
                      dest='outname',
                      action='store',
                      default=None,
                      help='The nwb filename')


    options,args = parser.parse_args()

    print(options)
    print(args)
    if ((options.subject is None)==(options.subject_id is None)):
        raise ValueError('Need to pass a subject with an ID')

    subject = get_subject_info(options.subject,options.subject_id)
    convert_NWB(f=args[0],
                exp_yaml=options.exp_yaml,
                subject=subject,
                NWBfilename=options.outname)
#
#TODO: use concatenated awake recordings to create an NWB file
#TODO: make the NWB file readable by KlustaKwik


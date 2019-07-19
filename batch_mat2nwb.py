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
                      help='A csv with the subject information')
    parser.add_option('-i','--id',
                      type='string',
                      dest='subject_id',
                      action='store',
                      default=None,
                      help='The ID number of the subject to grab from the csv')


    options,args = parser.parse_args()

    print(options)
    print(args)
    if ((options.subject is None)!=(options.subject_id is None)):
        raise ValueError('Need to pass a subject with an ID')

    if options.subject is None:
    # If no subject info is given, ask user for ID and genotype
        subject={}
        subject['ID'] = input('Enter the subject ID number')
        subject['genotype'] = input('Enter the subject genotype (default=c57bl6)')
        if len(subject['genotype'])==0:
            subject['genotype'] = 'c57bl6'
        subject = pynwb.file.Subject(subject_id=subject['ID'],genotype=subject['genotype'])
    else:
        subject = get_subject_info(options.subject,options.subject_id)

    file_list = glob.glob(os.path.join(args[0],'*.mat'))
    print(*file_list)
    print(subject)
    for f in file_list:
        convert_NWB(f=f,
                    exp_yaml=options.exp_yaml,
                    subject=subject)
#
#TODO: use concatenated awake recordings to create an NWB file
#TODO: make the NWB file readable by KlustaKwik


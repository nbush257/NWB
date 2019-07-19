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
    parser.add_option('-c','--concatenate',
                      dest='concat',
                      action='store_true',
                      default=False,
                      help='Flag to concatenate all the data files')

    parser.add_option('-i','--ignore_calib',
                      dest='ignore',
                      action='store_true',
                      default=False,
                      help='Flag to ignore calibration files')

    options,args = parser.parse_args()


    if options.subject is None:
    # If no subject info is given, ask user for ID and genotype
        print('\n')
        print('='*40)
        subject = input('Enter the subject ID number: ')

    file_list = glob.glob(os.path.join(args[0],'*.mat'))
    prefix = os.path.commonprefix(file_list)+'.nwb'
    if options.ignore:
        file_list = [i for i in file_list if 'calib' not in i]
    nwb_list = [os.path.splitext(i)[0]+'.nwb' for i in file_list]
    print('='*40)
    print(*file_list,sep='\n')
    print(subject)
    for f in file_list:
        convert_NWB(f=f,
                    exp_yaml=options.exp_yaml,
                    ID=subject)

    if options.concat:
        baselist = [os.path.split(i)[1] for i in file_list]
        concatenate_NWB(nwb_list,os.path.join(args[0],prefix))
        cleanup_nwb = input('Do you want to delete the temporary nwb files? ([Y]/n)')
        cleanup_mat = input('Do you want to delete the temporary matlab files? ([Y]/n)')

        if cleanup_mat.lower()=='y' or cleanup_mat=='':
            file_list = glob.glob(os.path.join(args[0],'*.mat'))
            for f in file_list:
                os.remove(f)
        if cleanup_nwb.lower()=='y' or cleanup_nwb=='':
            file_list = glob.glob(os.path.join(args[0],'*.mat'))
            for f in file_list:
                try:
                    nwb_f = os.path.splitext(f)[0]+'.nwb'
                    os.remove(nwb_f)
                except:
                    pass




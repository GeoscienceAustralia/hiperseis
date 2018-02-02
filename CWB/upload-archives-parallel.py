import boto3
import os
import multiprocessing
from botocore.utils import calculate_tree_hash
import math
import tarfile

base_location = '/g/data/ha3/GASeisDataArchive'

def ceiling_log(num, base):
    if num < 0:
       raise ValueError("Non-negative number only.")

    if num == 0:
       return 0

    return base ** (int(math.log(num, base))+1)

def transfer_part(transfer_file, start, sz, upload_id, vault_name):
    # got the error: PicklingError: Can't pickle <class 'botocore.client.Glacier'>: attribute lookup botocore.client.Glacier failed
    # when tried to pass the client as a paramter to each process per core
    # so have to create a local variable glacier_client
    glacier_client = boto3.client('glacier',
                                  region_name='ap-southeast-2',
                                  aws_access_key_id='DUMMY_ID',
                                  aws_secret_access_key='DUMMY_ACCESS')

    openfile = open(transfer_file, "rb")
    openfile.seek(start)
    data = openfile.read(sz)
    if (len(data) != sz):
        sz = len(data)
    partRange = "bytes {0}-{1}/*".format(start, (start+sz-1)) # Format '0-4194303'
    print("Sending Part:",partRange)
    ret = glacier_client.upload_multipart_part(vaultName=vault_name,
                                               uploadId=upload_id,
                                               range=partRange,
                                               body=data)
    if ret and 'ResponseMetadata' in ret and 'HTTPStatusCode' in ret['ResponseMetadata'] and ret['ResponseMetadata']['HTTPStatusCode'] == 204:
        return True
    else:
        return False

'''
def read_uploadid(file_name):
    with open('../upload_ids/'+file_name, 'r') as upload_id_infile:
        return upload_id_infile.read()

def write_uploadid(file_name, uploadId):
    with open('../upload_ids/'+file_name, 'w') as upload_id_outfile:
        upload_id_outfile.write(uploadId)
'''

def upload_archive(transfer_file, vault_name):
    glacier_client = boto3.client('glacier',
                                  region_name='ap-southeast-2',
                                  aws_access_key_id='DUMMY_ID',
                                  aws_secret_access_key='DUMMY_ACCESS')
    total = os.path.getsize(transfer_file)
    size = ceiling_log(total/32, 2) # 32 cores on this NCI machine

    init_mp_upl_resp = glacier_client.initiate_multipart_upload(vaultName=vault_name,
                                                                archiveDescription='2000_062 cwb waveform data',
                                                                partSize=str(size))
    print init_mp_upl_resp['uploadId']
    #write_uploadid(os.path.splitext(transfer_file)[0]+'.id', init_mp_upl_resp['uploadId'])

    '''
    # running 32 parallel thread with the code below caused a argument mismatch error in the
    # multiprocessing code! so continuing with serial upload for now. will investigate and
    # try to fix the parallel upload.
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    tasks = []
    start = 0
    while start < total:
        tasks.append( (transfer_file, start, size, init_mp_upl_resp['uploadId']) )
        start += size
    results = [pool.apply_async( transfer_part, t ) for t in tasks]

    all_uploaded=True
    for i, result in enumerate(results):
        if not result.get():
            print('Part number ', i, ' was not uploaded successfully')
            all_uploaded=False
    '''

    all_uploaded = True
    start = 0
    while start < total:
        if transfer_part(transfer_file, start, size, init_mp_upl_resp['uploadId'], vault_name):
            start += size
        else:
            all_uploaded = False

    if all_uploaded:
        print("All Files Uploaded")
        print("Verifying Checksum...")
        complete_up = glacier_client.complete_multipart_upload(vaultName=vault_name, uploadId=init_mp_upl_resp['uploadId'], archiveSize=str(total), checksum=calculate_tree_hash(open(transfer_file, 'rb')))
        print("Upload Completed:", complete_up)
    else:
        print("Upload of archive file:", transfer_file, " failed...")

def create_vault(year):
    glacier_client = boto3.client('glacier',
                                  region_name='ap-southeast-2',
                                  aws_access_key_id='DUMMY_ID',
                                  aws_secret_access_key='DUMMY_ACCESS')
    create_vault_resp = glacier_client.create_vault(vaultName=year)
    print create_vault_resp

def upload_yearly_archives(year):
    archives_location = os.path.join(base_location, year)
    if not os.path.isdir(archives_location):
        print('The folder for the year: ', year, ' does not exist at the location: ', base_location)
        return False

    files = filter(os.path.isfile, [os.path.join(archives_location, f) for f in os.listdir(archives_location)])
    no_of_files = len(files)
    files = [os.path.splitext(os.path.basename(f))[0] for f in files]
    files = list(set(files))
    no_of_days = len(files)
    # begin-section: below section is written to work around the abrupt killing of running process by NCI
    existing_files = [f for f in os.listdir('.') if f.endswith('.tar.gz')]
    existing_files = [os.path.splitext(os.path.splitext(f)[0])[0] for f in existing_files]
    files = [f for f in files if f not in existing_files]
    # end-section
    if no_of_files % 2 != 0 or no_of_days != no_of_files/2:
        print('There is something wrong with the expected waveform files\' number for the year: ', year)
        print('Please check that for each day, there are 2 files, one with .idx extension and another with .ms extension.')
        return False

    for file_name in files:
        tar = tarfile.open(file_name+'.tar.gz', "w:gz")
        tar.add(os.path.join(archives_location, file_name+'.idx'))
        tar.add(os.path.join(archives_location, file_name+'.ms'))
        tar.close()
        upload_archive(file_name+'.tar.gz', 'cwb-archive-'+year)

#upload_archive('2000_062_d.tar.gz')
create_vault('cwb-archive-2004')
upload_yearly_archives('2004')


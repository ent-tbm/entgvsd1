import zlib
import sys
import os
import subprocess

def count_gz(fname):
    """Counts the number of .gz extensions on a filename"""
    pieces = fname.split('.')
    for i in range(len(pieces)-1,-1,-1):
        if pieces[i] != 'gz' or i==0:
            return '.'.join(pieces[:i+1]), len(pieces)-1-i    

def copy_recompress(iroot, oroot, dir, zleaf):
    leaf,ngz = count_gz(zleaf)

    # Don't copy some files
    if leaf == '.DS_Store':
        return

    oname = os.path.join(oroot, dir, leaf)
    if os.path.exists(oname):
        # Assume the file was created properly (atomically)
        # So skip it
        return

    print('Writing {}'.format(oname))
    sys.stdout.flush()

    # If file doesn't need decompressing, just symlink it
    if ngz == 0:
        os.symlink(os.path.join(iroot,dir,zleaf), oname)
        return

    # Decompress the .gz file
    step = 0
    dcs = []
    for x in range(0,ngz):
        dcs.append(zlib.decompressobj(32 + zlib.MAX_WBITS))

    fname = os.path.join(iroot,dir,zleaf)
    with open(fname, 'rb') as zin:
        with open(os.path.join(oroot, '_TMP_{}'.format(step)), 'wb') as fout:
            step += 1
            while True:
                data = zin.read(1024*1024)   # Too big causes segfault
                if len(data) == 0:
                    break

                # Decompress chunk in stages
                for i in range(0,ngz):
                    data1 = dcs[i].decompress(data)
                    data = data1

                fout.write(data)

            # https://stackoverflow.com/questions/2333872/atomic-writing-to-file-with-python
            fout.flush()
            os.fsync(fout.fileno())

    # Recompress if it's NetCDF
    if leaf.endswith('.nc'):
        cmd = ['nccopy', '-k', 'nc4', '-d', '4', '-s',
            os.path.join(oroot, '_TMP_{}'.format(step-1)),
            os.path.join(oroot, '_TMP_{}'.format(step))]
        step += 1
        err = subprocess.call(cmd)
        if (err != 0):
            sys.stderr.write('Error {}\n'.format(os.path.join(iroot,dir,zleaf)))
            sys.stderr.flush()
            return

    # Move to final location
    os.rename(
        os.path.join(oroot, '_TMP_{}'.format(step-1)),
        oname)

def copy_all(iroot, oroot):

    n=0
    for idir, subdirList, fileList in os.walk(iroot):
        dir = idir[len(iroot)+1:]
        os.makedirs(os.path.join(oroot,dir), exist_ok=True)

        for fname in fileList:
            copy_recompress(iroot, oroot, dir, fname)
            n += 1


#IROOT = '/home2/rpfische/entgvsd0_orig'
OROOT = '/home2/rpfische/entgvsd0'
IROOT = '/home2/rpfische/entgvsd0_orig'

copy_all(IROOT,OROOT)


sys.exit(0)





# https://stackoverflow.com/questions/12571913/python-unzipping-stream-of-bytes
dc0 = zlib.decompressobj(32 + zlib.MAX_WBITS)
dc1 = zlib.decompressobj(32 + zlib.MAX_WBITS)

with open('msg.gz.gz') as zzin:
    # https://stackoverflow.com/questions/519633/lazy-method-for-reading-big-file-in-python
    while True:
        datazz = zzin.read(1024*10)
        if len(datazz) == 0:
            break

        dataz = dc0.decompress(datazz)
        data = dc1.decompress(dataz)
        sys.stdout.write(data)
        #print(len(datazz),len(dataz))

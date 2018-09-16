#!/bin/env python3
#

import netCDF4
import os
import sys
import subprocess
import zlib
import time
import traceback

# Size of our Chunk grid
IM_CHUNK=18
JM_CHUNK=15

def decompress_gzz(iname, ngz, oname):

    pipeline = []
    pipeline.append(('zcat', iname))
    for i in range(0,ngz-1):
        pipeline.append(('zcat'))

    with open(oname, 'wb') as out:
        procs = []
        procs.append(subprocess.Popen(pipeline[0], stdout=subprocess.PIPE))
        for cmd in pipeline[1:-1]:
            procs.append(subprocess.Popen(cmd, stdin=procs[-1].stdout, stdout=subprocess.PIPE))
        procs.append(subprocess.Popen(pipeline[-1], stdin=procs[-1].stdout, stdout=out))
        procs[-1].communicate()


def copy_recompress(iroot, oroot, dir, leaf):

    print('=============================')
    print('iroot',iroot)
    print('oroot',oroot)
    print('dir',dir)
    print('leaf',leaf)

    os.makedirs(os.path.join(oroot,dir), exist_ok=True)

    # Find the least-compressed version of this file
    iname = os.path.join(iroot, dir, leaf)
    ngz = 0
    while True:
        if os.path.exists(iname):
            break
        iname += '.gz'
        ngz += 1
        if ngz > 2:
            raise IOError('Cannot find input file {}'.format(os.path.join(iroot,dir,leaf)))

    oname = os.path.join(oroot, dir, leaf)
    if os.path.exists(oname):
        # Assume the file was created properly (atomically)
        # So skip it
        return

    print('Writing {}'.format(oname))
    sys.stdout.flush()

    # If file doesn't need decompressing, just symlink it
    if ngz == 0:
        os.symlink(iname, oname)
        return


    step=0
    decompress_gzz(iname, ngz, os.path.join(oroot, '_TMP_{}'.format(step)))
    step += 1


#    # Decompress the .gz file
#    step = 0
#    dcs = []
#    for x in range(0,ngz):
#        dcs.append(zlib.decompressobj(32 + zlib.MAX_WBITS))
#
#    zlib_chunk_size = 8*1024   # Too large doesn't work
#    sizemb = os.path.getsize(iname) / float(zlib_chunk_size)
#    sys.stdout.write('    ...unzipping ({:0.1f} Mb)'.format(sizemb))
#    sys.stdout.flush()
#    t0 = time.time()
#    with open(iname, 'rb') as zin:
#        with open(os.path.join(oroot, '_TMP_{}'.format(step)), 'wb') as fout:
#            step += 1
#            while True:
#                data = zin.read(zlib_chunk_size)   # Too big causes segfault
#                if len(data) == 0:
#                    break
#
#                # Decompress chunk in stages
#                for i in range(0,ngz):
#                    data1 = dcs[i].decompress(data)
#                    data = data1
#
#                fout.write(data)
#                sys.stdout.write('.')
#                sys.stdout.flush()
#
#            # https://stackoverflow.com/questions/2333872/atomic-writing-to-file-with-python
#            fout.flush()
#            os.fsync(fout.fileno())
#
#            t1 = time.time()
#            print(' [{:0.0f}s]'.format(t1-t0))
#



# Don't recompress.  It's slow, and nccopy doesn't really work
# Need to make a Python program that recompresses.
# No real need to recompress anyway
#    # ---------- Recompress if it's NetCDF
#    if leaf.endswith('.nc'):
#        with netCDF4.Dataset(os.path.join(oroot, '_TMP_{}'.format(step-1))) as nc:
#            chunk_specs = []
#            chunk_specs.append('lon/{:0.0f}'.format(len(nc.dimensions['lon']) / IM_CHUNK))
#            chunk_specs.append('lat/{:0.0f}'.format(len(nc.dimensions['lat']) / JM_CHUNK))
#            done_dims = {'lon', 'lat'}
#            for dname,_ in nc.dimensions.items():
#                if dname not in done_dims:
#                    chunk_specs.append('{}/1'.format(dname))
#
#        sys.stdout.write('   ...converting to chunked NetCDF4 {}'.format(','.join(chunk_specs)))
#        sys.stdout.flush()
#        cmd = ['nccopy', '-w', '-k', 'nc4', '-d', '1', '-s',
#            '-c', ','.join(chunk_specs),
#            os.path.join(oroot, '_TMP_{}'.format(step-1)),
#            os.path.join(oroot, '_TMP_{}'.format(step))]
#        step += 1
#        err = subprocess.call(cmd)
#        t2 = time.time()
#        print(' [{:0.0f}s]'.format(t1-t0))
#        if (err != 0):
#            sys.stderr.write('Error {}\n'.format(os.path.join(iroot,dir,leaf)))
#            sys.stderr.flush()
#            return

    # Move to final location
    os.rename(
        os.path.join(oroot, '_TMP_{}'.format(step-1)),
        oname)



IROOT=sys.argv[1]
OROOT=sys.argv[2]
DIR=sys.argv[3]          # Directory of the file (after the root)
LEAF=sys.argv[4]         # Leaf name (without the .nc)

copy_recompress(IROOT, OROOT, DIR, LEAF)

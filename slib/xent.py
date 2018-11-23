#!/usr/bin/env python
#

from __future__ import print_function
import os
import sys
import subprocess

# --------------------------------------------------
def search_up(path, condition_fn):
    """Scans up a directory tree until a condition is found on one of the paths."""

    path = os.path.abspath(path)#os.path.split(fname)[0])
    while True:
        ret = condition_fn(path)
        if ret is not None:
            return ret
        new_path = os.path.dirname(path)
        if new_path == path:
            # Didn't find it.
            return None
        path = new_path

def has_file(path, fname):
    ret = os.path.join(path, fname)
    if os.path.exists(ret):
        return ret
    return None
# --------------------------------------------------

# Scan for first non-flag argument, that will be script (Fortran source) name
for i in range(1,len(sys.argv)):
    if len(sys.argv[i]) > 0 and sys.argv[i][0] != '-':
        fsrc = os.path.abspath(sys.argv[i])

        # Find EntGVSD root
        root = os.path.split(search_up(os.path.split(fsrc)[0], lambda x : has_file(x, 'entgvsd_root.txt')))[0]
        if root is None:
            raise ValueError('Source file {} is not in an EntGVSD source directory!'.format(fsrc))

        # Use build/ directory
        launcher = os.path.join(root, 'build', 'bin', 'entgvsd')

        cmd = [launcher] + sys.argv[1:]
        print(cmd)
        sys.exit(subprocess.call(cmd))

raise ValueError('Could not identify Fortran source file on command line')

import sys
import os
import subprocess

# When run from a data/ directory (or subdirectory)... recursively
# downloads data that belongs there.

this_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path = [this_dir] + sys.path

import nasaportal

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

def any_isfile(fnames):
    """Returns the filename if any exist; otherwise none"""
    for fname in fnames:
        if os.path.isfile(fname):
            return fname
    return None


def _condition(path):
    root,leaf = os.path.split(path)
    if leaf=='data':
        return root
    else:
        return None

# ------ Get directory on NASA portal corresponding to this one
# Search up till we find 'data'
cwd = os.path.abspath('.')
project_root = search_up(cwd, _condition)
portal_dir = '/GISS/Ent_TBM/EntGVSD' + cwd[len(project_root):]

sys.stderr.write('Looking for files to download...\n')
for root, dirs, files in nasaportal.walk(portal_dir):
    for name in files:
        iname = '/'.join((root, name))
        oname = '.' + iname[len(portal_dir):]

        # Possible output names to look for
        possible_onames = [oname]
        if oname.endswith('.gz'):
            possible_onames.append(oname[:-3])

        exists = any_isfile(possible_onames)
        if (exists is not None):
            sys.stderr.write('File exists: {}\n'.format(exists))
        else:
            sys.stderr.write('Downloading {}\n'.format(oname))
            nasaportal.download(iname, ofile=oname)

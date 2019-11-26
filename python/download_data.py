import sys
import os
import subprocess
import argparse

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

def do_main(action):
    """clean: bool
        If set, we want to DELETE downloaded files, instead of downloading them.
    """
    # ------ Get directory on NASA portal corresponding to this one
    # Search up till we find 'data'
    cwd = os.path.abspath('.')
    project_root = search_up(cwd, _condition)
    portal_dir = '/GISS/Ent_TBM/EntGVSD' + cwd[len(project_root):]

    sys.stderr.write('Looking for files to {}...\n'.format(action))
    for root, dirs, files in nasaportal.walk(portal_dir):
        for name in files:
            iname = '/'.join((root, name))
            oname = '.' + iname[len(portal_dir):]

            # Possible output names to look for
            possible_onames = [oname]
            if oname.endswith('.gz'):
                possible_onames.append(oname[:-3])

            if action == 'list':
                exists = any_isfile(possible_onames)
                if exists is not None:
                    print('* {}'.format(exists))
                else:
                    print('  {}'.format(oname))
            elif action == 'clean':
                # We want to delete downloaded files, not download them
                for on in possible_onames:
                    try:
                        os.remove(on)
                        print('Removed {}'.format(on))
                    except FileNotFoundError:
                        pass
            elif action == 'download':
                exists = any_isfile(possible_onames)
                if (exists is not None):
                    sys.stderr.write('File exists: {}\n'.format(exists))
                else:
                    sys.stderr.write('Downloading {}\n'.format(oname))
                    nasaportal.download(iname, ofile=oname)
            else:
                raise ValueError('Illegal action: {}'.format(action))

def main():
    parser = argparse.ArgumentParser(description='Download corresponding EntGVSD data files')
    parser.add_argument('action', type=str, default='download',
        help='Action to take: list, download, or clean')
#    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
#        help='DELETE data files instead of downloading them')
    args = parser.parse_args()
    do_main(args.action)


main()

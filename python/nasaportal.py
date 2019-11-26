#import requests
import subprocess
import re
import os

# https://www.nccs.nasa.gov/nccs-users/instructional/ftp-deprecation

hrefRE = re.compile(r'<a\s+href="([^"]*)">', flags=re.MULTILINE)

with open('/home2/rpfische/tmp/y') as fin:
    str = fin.read()


PORTAL_ROOT = 'https://portal.nccs.nasa.gov/datashare/'
def list_hrefs(dir):
    url = PORTAL_ROOT + dir
    if url[-1] != '/':
        url = url + '/'
    cmd = ['curl', '--location', url]
    with open(os.devnull, 'w') as FNULL:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=FNULL)
    html = proc.stdout.read().decode()

    for matches in hrefRE.finditer(html):
#        requests.get(url, allow_redirects=True).content):
        yield matches.group(1)

def list_files(dir):
    state = 0
    for href in list_hrefs(dir):
        if state==0:
            if href.startswith('/datashare'):
                state=1
                continue
        if state==1:
            if href.startswith('http://'):
                return
            yield href

def walk(dir):
    es


print(list(list_files('GISS/Ent_TBM/EntGVSD/inputs/lai/BNU/doy/2004')))

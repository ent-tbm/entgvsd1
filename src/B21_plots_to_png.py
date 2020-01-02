import subprocess
import os

"""Convert R-generated PDF plots to PNG format, for easy display"""

NTASKS=10    # Number of Ghostscript tasks to spawn at once
root = '../outputs'
idir = os.path.join(root, 'plots')
odir = os.path.join(root, 'png')

os.makedirs(odir, exist_ok=True)

ileaves = sorted(os.listdir(idir))

print('---------------- Already Processed')
jobs = []
for ileaf in ileaves:
    # Ignore extraneous files
    if not ileaf.endswith('.pdf'):
        continue

    ifname = os.path.join(idir, ileaf)
    oleaf = ileaf
    oleaf = oleaf.replace('_forplot.nc.pdf', '.png')
    oleaf = oleaf.replace('.nc.pdf', '.png')
    oleaf = oleaf.replace('.pdf', '.png')
    ofname = os.path.join(odir, oleaf)

    itime = os.path.getmtime(ifname)
    otime = os.path.getmtime(ofname) if os.path.exists(ofname) else 0

    if itime > otime:
        jobs.append((ifname, ofname))
    else:
        print(ofname)

processes = []
for i in range(0,len(jobs),NTASKS):
    xjobs = jobs[i:min(len(jobs),i+NTASKS)]

    # Spawn NTASKS processes
    print('--------------- Converting {}-{} (of {})'.format(i+1,i+len(xjobs), len(jobs)))
    for ifname,ofname in xjobs:
        print(ofname)
        cmd = ['gs', '-dQUIET', '-dUseCropBox', '-dBATCH', '-dNOPAUSE', '-r600', '-sPAPERSIZE=letter', '-sDEVICE=png16m', '-sOutputFile='+ofname, ifname]
        processes.append(subprocess.Popen(cmd))

    # Wait for all NTASKS process to finish
    for process in processes:
        process.wait()

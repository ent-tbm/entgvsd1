import subprocess

cmd = ('python3', 'plotregion.py')

while True:
    err = subprocess.call(cmd)
    if err == 0:
        break

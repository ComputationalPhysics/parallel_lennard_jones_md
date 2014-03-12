import subprocess

subprocess.call('rm *.txt', shell = True)
subprocess.call('rm -rf state_files', shell = True)
subprocess.call('rm -rf movie_files', shell = True)
subprocess.call('rm md.ini', shell = True)
subprocess.call('rm *.o', shell = True)
subprocess.call('rm log', shell = True)
subprocess.call('rm *.pyc', shell = True)
subprocess.call('rm Tocontinue', shell = True)

print "Cleanup complete"
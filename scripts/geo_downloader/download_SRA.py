# Make download links for GEO, querying ncbi ftp
# Zijun Zhang, 10.26.2016

import sys
import os
from ftplib import FTP
import yaml

ftp = FTP('ftp-trace.ncbi.nlm.nih.gov')
ftp.login()

#dir = sys.argv[1]
project = sys.argv[1]
#dir = os.path.abspath(dir)
dir = os.path.join('projects', project, 'sra')
with open(os.path.join('projects', project, 'config', 'config.yaml'), 'r') as f:
	config = yaml.load(f)

sample_type_dict = config['sample_type_dict']

for line in sys.stdin:
	ft, ftp_dir_path = line.rstrip().split()
	ftp.cwd('/')
	srx_dir = '/'.join(ftp_dir_path.split("/")[3:])
	ftp.cwd(srx_dir)
	sra_dir = []
	ftp.dir(sra_dir.append)
	sra = [x.split()[-1] for x in sra_dir]
	try:
		type = sample_type_dict[ft]
	except KeyError:
		continue
	for s in sra:
		print >> sys.stdout, 'wget -O {0} {1}'.format(dir+'/'+ type+'/' + ft+'.sra', '/'.join([ftp_dir_path, s, s + '.sra']))
		
	

ftp.quit()

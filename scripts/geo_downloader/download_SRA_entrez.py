## use Bio::Entrez to download the SRR 
## to a structured folder
## Zijun Zhang
## 3.30.2018

import sys
import os
import json
import time
import pandas as pd
import yaml
from Bio import Entrez
from collections import defaultdict
from GEOparse import utils

try:
    from urllib.error import HTTPError
except ImportError:
    from urllib2 import HTTPError

EMAIL = "admin@local.host"

Entrez.email = EMAIL

sample_to_url = defaultdict(list)

for line in sys.stdin:
	ele = line.strip().split()
	sample_name, srx = ele[0], ele[1]
	searchdata = Entrez.esearch(db='sra', term=srx, usehistory='y',
								retmode='json')
	answer = json.loads(searchdata.read())
	ids = answer["esearchresult"]["idlist"]
	assert len(ids) == 1, "There should be one and only one ID per SRX"
	
	# using ID fetch the info
	number_of_trials = 10
	wait_time = 30
	for trial in range(number_of_trials):
		try:
			results = Entrez.efetch(db="sra", id=ids[0],
									rettype="runinfo",
									retmode="text").read()
			break
		except HTTPError as httperr:
			if "502" in str(httperr):
				print(("Error: %s, trial %i out of %i, waiting "
							  "for %i seconds.") % (
								 str(httperr),
								 trial,
								 number_of_trials,
								 wait_time))
				time.sleep(wait_time)
			else:
				raise httperr
	
	df = pd.DataFrame(
		[i.split(',') for i in results.split('\n') if i != ''][1:],
		columns=
		[i.split(',') for i in results.split('\n') if i != ''][0])
	
	sample_to_url[sample_name] = [x for x in df['download_path'].values ]


project = sys.argv[1]
dir = os.path.join('projects', project, 'sra')
with open(os.path.join('projects', project, 'config', 'config.yaml'), 'r') as f:
	config = yaml.load(f)

sample_type_dict = config['sample_type_dict']

for sample_name in sample_to_url:
	if not sample_name in config['sample_type_dict']:
		continue
	type = config['sample_type_dict'][sample_name]
	for s in sample_to_url[sample_name]:
		ft = os.path.basename(s)
		#print >> sys.stdout, 'mkdir -p {0}'.format(
		#	os.path.join(dir, type, sample_name))
		filepath = os.path.join(dir, sample_name, ft+'.sra')
		directory_path = os.path.join(dir, sample_name)
		print >> sys.stdout, 'wget -O {0} {1}'.format(
			filepath, 
			s)
		utils.mkdir_p(os.path.abspath(directory_path))
		url = s
		utils.download_from_url(url, filepath, aspera=False)
		
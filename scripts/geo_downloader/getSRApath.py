## the Python equivalent of getSRA
## Zijun Zhang
## 3.29.2018

import sys
import re
import GEOparse
from collections import defaultdict

gid = sys.argv[1]
outfp = sys.argv[2]
gse = GEOparse.get_GEO(gid, destdir="/tmp")

queries = defaultdict(list)
for gsm_name in gse.gsms:
	gsm = gse.gsms[gsm_name]
	title = gsm.metadata['title'][0]
	title = re.sub(r"[\s|\-|/]","_", title)
	title = re.sub(r"[,|\+]","", title)
	for sra in gsm.relations['SRA']:
		query = sra.split("=")[-1]
		assert 'SRX' in query, "Sample looks like it is not SRA: %s" % query
		queries[title].append(query)

ftpaddres = ("ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/"
                     "ByExp/sra/SRX/{range_subdir}/{file_dir}/"
                     "{file_dir}")

with open(outfp, 'w') as f:
	for sample_name in queries:
		assert len(queries[sample_name])==1, "More than one SRA associated to %s"% sample_name
		srx = queries[sample_name][0]
		url = ftpaddres.format(range_subdir=srx[0:6], file_dir=srx)
		#f.write("%s\t%s\n"%(sample_name, url))
		f.write("%s\t%s\n"%(sample_name, srx))

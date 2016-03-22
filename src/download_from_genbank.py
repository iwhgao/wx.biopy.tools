#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Email: gaowenhui2012@gmail.com
# Version: v1.0.0
# Created: 2016/3/22
# ---------------------------


import sys
import os
import urllib2
from Bio import SeqIO
from Bio import Entrez
from wxbiopyutils import parse_config, gen_logger


def download_from_genbank(id_file, output_file, db="nucleotide", rettype="fasta", retmode="text", save_format="fasta"):
	"""
	:type id_file: basestring
	:type output_file: basestring
	"""

	Entrez.email = 'A.N.Other@example.com'
	id_file_handle = open(id_file, 'r')
	id_lines = id_file_handle.readlines()
	id_lines = [id.strip('\t') for id in id_lines]
	id_list_str = ','.join(id_lines)

	handle = Entrez.efetch(db=db, rettype=rettype, retmode=retmode, id=id_list_str)

	output_file_handle = open(output_file, 'w')
	cnt = 1
	for seq_record in SeqIO.parse(handle, rettype):
		SeqIO.write(seq_record, output_file_handle, save_format)
		log.info('#{1} Processed {0}'.format(seq_record.id, cnt))
		cnt += 1
	handle.close()
	log.info('Total {0} queries done!'.format(cnt - 1))


if __name__ == '__main__':
	log_basename = (os.path.split(sys.argv[0]))[1]
	log = gen_logger(log_basename)
	log.info('Start!')

	conf = parse_config('../conf/download_from_genbank.ini')

	# 代理
	proxy = 'web-proxy.oa.com:8080'
	opener = urllib2.build_opener(urllib2.ProxyHandler({'http': proxy}))
	urllib2.install_opener(opener)

	if conf and conf['idlist']['file_path'] and conf['output']['file_path'] and conf['output']['format']:
		download_from_genbank(conf['idlist']['file_path'],
							  conf['output']['file_path'],
							  db=conf['genbank']['db'],
							  rettype=conf['genbank']['rettype'],
							  retmode=conf['genbank']['retmode'],
							  save_format=conf['output']['format'])
		log.info('Read configure file')
	else:
		log.error('Read configure file error')

	log.info('Finished!')

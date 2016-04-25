#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Email: gaowenhui2012@gmail.com
# Version: v1.0.0
# Created: 2016/3/22
# ---------------------------


import urllib2
from Bio import ExPASy
from Bio import SeqIO
from wxbiopyutils import parse_config, gen_logger, basedir


def download_from_swissprot(id_file, output_file, rettype="swiss", save_format="swiss"):
	"""
	:type id_file: basestring
	:type output_file: basestring
	"""

	cnt = 1
	output_file_handle = open(output_file, 'w')
	with open(id_file, 'r') as f:
		for line in f:
			query_id = line.strip()

			try:
				handle = ExPASy.get_sprot_raw(query_id)
			except urllib2.HTTPError as e:
				log.warning('{0} query failed'.format(query_id))

			seq_record = SeqIO.read(handle, rettype)
			SeqIO.write(seq_record, output_file_handle, save_format)
			log.info('#{1} Processed {0}'.format(seq_record.id, cnt))
			cnt += 1
			handle.close()

	log.info('Total {0} queries done!'.format(cnt - 1))

if __name__ == '__main__':
	log_basename = __file__
	log = gen_logger(log_basename)
	log.info('Start!')

	conf = parse_config(basedir + '/conf/download_from_swissprot.ini')

	# 代理
	if 'is_proxy' in conf['global'] and conf['global']['is_proxy'] == "1":
		proxy = conf['proxy']['proxy']
		opener = urllib2.build_opener(urllib2.ProxyHandler({'http': proxy}))
		urllib2.install_opener(opener)

	if conf \
			and conf['idlist']['file_path'] \
			and conf['output']['file_path'] \
			and conf['swissprot']['rettype'] \
			and conf['output']['format']:
		download_from_swissprot(conf['idlist']['file_path'],
								conf['output']['file_path'],
								rettype=conf['swissprot']['rettype'],
								save_format=conf['output']['format'])
		log.info('Read configure file')
	else:
		log.error('Read configure file error')

	log.info('Finished!')


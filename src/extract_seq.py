#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Version: v1.0.0
# Created: 2016/3/21
# ---------------------------

import sys
import os
import re
from Bio import SeqIO
from Bio import SeqRecord
from wxbiopyutils import parse_config, gen_logger


def extract_seq_from_file(seq_file, coords_file, output_file):
	# 记录reference sequence名称
	chrs = []

	# 存储片段
	chr_seg = {}

	# 对片段计数
	cnt = 0

	seqio = SeqIO.parse(seq_file, 'fasta')
	for seq_record in seqio:
		chrs.append(seq_record.id)

	with open(coords_file, 'r') as f:
		for line in f:
			cnt += 1
			line = line.strip('\n')
			regions = re.split('\s+', line)

			if regions[0] not in chrs:
				log.warning('{0} not in reference sequence'.format(regions[0]))

			if len(regions) < 3:
				log.warning('The numbers of this line are less than 3(required)')
				continue

			if regions[0] not in chr_seg:
				chr_seg[regions[0]] = []
				chr_seg[regions[0]].append(regions)
			else:
				chr_seg[regions[0]].append(regions)

	log.info('Summary: {0} chromosomes, {1} segments processed'.format(len(chr_seg), cnt))

	res_file_handle = open(output_file, 'w')

	# 遍历reference sequence
	seqio = SeqIO.parse(seq_file, 'fasta')
	for seq_record in seqio:
		if seq_record.id in chr_seg:
			for seg in chr_seg[seq_record.id]:
				try:
					# 创建SeqRecord对象
					tmp_seq = SeqRecord.SeqRecord(seq=(seq_record.seq)[(int(seg[1])-1):int(seg[2])],
												  id='{0}:{1}..{2}:{3}'.format(seg[0], seg[1], seg[2], seg[3]))

					# 当strang为-时， 进行反向互补处理
					if seg[3] == '-':
						tmp_seq = tmp_seq.reverse_complement(id=True,
															 name=True,
															 description='reverse_complement')

					SeqIO.write(tmp_seq, res_file_handle, 'fasta')
				except Exception as e:
					log.error(e)
		else:
			log.warning(seq_record.id + ' not exists in reference sequences')

	res_file_handle.close()


if __name__ == '__main__':

	log_basename = (os.path.split(sys.argv[0]))[1]
	log = gen_logger(log_basename)
	log.info('Start!')

	conf = parse_config('../conf/extract_seq.ini')

	if conf and conf['seq']['file_path'] and conf['coords']['file_path'] and conf['output']['file_path']:
		extract_seq_from_file(conf['seq']['file_path'], conf['coords']['file_path'], conf['output']['file_path'])
	else:
		log.error('Configure file error!')

	log.info('Finished!')

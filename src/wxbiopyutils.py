#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Email: gaowenhui2012@gmail.com
# Version: v1.0.0
# Created: 2016/3/21
# ---------------------------


import os
import re
import logging
from configobj import ConfigObj


# 获取src的上一级目录
basedir = os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def parse_config(conf_file):
	"""
	解析配置文件
	:param conf_file:
	:return:
	"""
	config = None
	try:
		config = ConfigObj(conf_file, encoding='UTF8')
	except Exception as e:
		print e

	return config


def gen_logger(log_file):
	"""
	获取log对象
	:param log_file:
	:return:
	"""
	log_file = re.sub('\.py$', '.log', log_file)
	log_file = re.sub('\.exe$', '.log', log_file)

	log_file = basedir + '/log/' + log_file

	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s|%(levelname)s|%(message)s',
						datefmt='%Y-%m-%d %H:%M:%S',
						filename=log_file,
						filemode='a')
	log = logging.getLogger()

	return log


if __name__ == '__main__':
	print basedir

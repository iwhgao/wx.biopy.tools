#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Email: gaowenhui2012@gmail.com
# Version: v1.0.0
# Created: 2016/3/21
# ---------------------------


from configobj import ConfigObj
import logging

def parse_config(conf_file):

    config = None
    try:
        config = ConfigObj(conf_file, encoding='UTF8')
    except Exception as e:
        print e

    return config


def gen_logger(log_file):

    log_file = log_file.replace('.py', '.log')
    log_file = log_file.replace('.exe', '.log')

    log_file = '../log/' + log_file

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s|%(levelname)s|%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=log_file,
                        filemode='a')
    log = logging.getLogger()

    return log

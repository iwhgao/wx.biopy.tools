#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Email: gaowenhui2012@gmail.com
# Version: v1.0.0
# Created: 2016/3/10
# ---------------------------


from distutils.core import setup
import py2exe


py2exe_options = {
	"dll_excludes": ["w9xpopen.exe", "numpy-atlas.dll"],
	"compressed": 1,
	"optimize": 2,
	"ascii": 0,
	"bundle_files": 1,
}

setup(
	name='Bioinformatics python tools on windows platform',
	description='Bioinformatics python tools on windows platform',
	version='1.0.0',
	console=[{'script': 'drawChromosome4Cotton.py'}],
	zipfile=None,
	options={'py2exe': py2exe_options},
	author='deangao',
	author_email='gaowenhui2012@gmail.com'
)

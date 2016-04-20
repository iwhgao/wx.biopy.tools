#!/usr/bin/python
# -*- coding:utf-8 -*-

# ---------------------------
# Author: deangao
# Copyright: 2016 deangao
# Version: v1.0.0
# Created: 2016/3/14
# ---------------------------


import re
import sys
import logging
import svgwrite
from configobj import ConfigObj


def parse_config(conf_file):
	"""
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
	:param log_file:
	:return:
	"""

	log_file = log_file.replace('.py', '.log')
	log_file = log_file.replace('.exe', '.log')

	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s|%(levelname)s|%(message)s',
						datefmt='%Y-%m-%d %H:%M:%S',
						filename=log_file,
						filemode='a')
	log = logging.getLogger()
	return log


def get_chr_layout(conffile):
	"""
	获取chr的层次及每层的顺序
	:param conffile:
	:return: layouts
	"""

	# TODO: 暂时默认弄两层， 后面再优化
	layouts = []
	layouts.append([])
	layouts.append([])
	chr_order_layout = conffile['chromosome_settings']['orders']
	for chr in chr_order_layout:
		res = re.findall('(.*?)\(([0-9]+)\)', chr)
		if len(res[0]) == 2:
			layouts[int(res[0][1]) - 1].append(res[0][0])

	# print layouts
	return layouts


def parse_chr_profile(chr_profile):
	"""
	# 获取chr的ID， name， length，color属性
	# 获取chrs中的最大长度
	:param chr_profile:
	:return: max_length
	:return chrs：w
	"""

	max_length = 0
	chrs = {}
	with open(chr_profile, 'r') as f:
		for line in f:
			if re.match('^#', line):
				continue

			line = re.sub('\n$', '', line)
			regions = re.split('\t|\s+', line)

			if not regions or len(regions) < 4:
				log.error('Number of fields in chromosome profile is less then 4')
			if max_length < int(regions[2]):
				max_length = int(regions[2])

			chrs[regions[0]] = regions

	return max_length, chrs


def get_chr_posi(conffile):
	layouts = get_chr_layout(conffile)
	max_length, chrs = parse_chr_profile(conffile['chromosome_settings']['chr_profile'])

	left_margin = int(conffile['chromosome_settings']['left_margin'])
	chr_left_margin = int(conffile['chromosome_settings']['chr_left_margin'])
	layout_top_margin = int(conffile['chromosome_settings']['layout_top_margin'])
	chr_label_height = int(conffile['chromosome_settings']['chr_label_height'])
	chr_height = int(conffile['chromosome_settings']['chr_height'])
	chr_width = int(conffile['chromosome_settings']['chr_width'])

	chr_posi = {}
	factor = float(chr_height) / max_length
	for l, lay_chrs in enumerate(layouts):
		for c, chr in enumerate(lay_chrs):
			chr_profile = {}
			chr_profile["x"] = left_margin + c * chr_width + c * chr_left_margin
			chr_profile["y"] = (l+1)*layout_top_margin + (l+1)*chr_label_height + l * chr_height
			chr_profile["w"] = chr_width
			chr_profile["h"] = int(float(chrs[chr][2]) * factor)
			chr_profile["chr"] = chrs[chr]
			chr_posi[chr] = chr_profile

	# print chr_posi
	return chr_posi, factor


def get_chr_label_posi(conffile):
	chr_label = {}
	all_chr_posi, fac = get_chr_posi(conffile)

	for chr_tmp in all_chr_posi:
		chr_label_elem = {
			"x": int(all_chr_posi[chr_tmp]["x"] + float(conffile['chromosome_settings']['chr_width']) / 2),
			"y": int(all_chr_posi[chr_tmp]["y"] - float(conffile['chromosome_settings']['chr_label_height']) / 2),
			"title": all_chr_posi[chr_tmp]["chr"][1]}
		chr_label[chr_tmp] = chr_label_elem

	return chr_label


def get_mark_label_posi(all_chrs, fac, mark_profile, conffile):
	marks_ticks = []
	marks_labels = []
	with open(mark_profile, 'r') as f:
		for line in f:
			mark_tick = {}
			mark_label = {}
			if re.match('^#', line):
				continue

			line = re.sub('\n$', '', line)
			regions = re.split(',', line)

			if not regions or len(regions) != 8:
				log.error('Number of fields in chromosome profile is less then 4')
			elif len(regions) == 8:
				if regions[7] == "" or regions[7] == "right":
					mark_tick["x1"] = int(all_chrs[regions[0]]["x"] + float(conffile['chromosome_settings']['chr_width']))
					mark_tick["y1"] = int(all_chrs[regions[0]]["y"] + float(regions[1])*fac)
					mark_tick["x2"] = mark_tick["x1"] + int(regions[4])
					mark_tick["y2"] = mark_tick["y1"]
					mark_tick["color"] = regions[6]
				elif regions[7] == "left":
					mark_tick["x1"] = int(all_chrs[regions[0]]["x"])
					mark_tick["y1"] = int(all_chrs[regions[0]]["y"] + float(regions[1])*fac)
					mark_tick["x2"] = mark_tick["x1"] - int(regions[4])
					mark_tick["y2"] = mark_tick["y1"]
					mark_tick["color"] = regions[6]

				marks_ticks.append(mark_tick)

				if regions[2] != "":
					mark_label["x"] = mark_tick["x2"]
					mark_label["y"] = mark_tick["y2"]
					mark_label["title"] = regions[2]
					mark_label["color"] = regions[5]
					marks_labels.append(mark_label)

	return marks_ticks, marks_labels


def draw_pic(dwg, conffile):
	rx = int(conf['chromosome_settings']['rx'])
	ry = int(conf['chromosome_settings']['ry'])
	all_chr_posi, fac = get_chr_posi(conffile)
	all_chr_label = get_chr_label_posi(conffile)
	all_marks_ticks, all_marks_labels = get_mark_label_posi(all_chr_posi, fac, conffile['chromosome_settings']['label_profile'], conffile)

	# 绘制chr框架
	for chr_tmp in all_chr_posi:
		dwg.add(svgwrite.shapes.Rect(insert=(all_chr_posi[chr_tmp]['x'], all_chr_posi[chr_tmp]['y']),
									 size=(all_chr_posi[chr_tmp]['w'], all_chr_posi[chr_tmp]['h']),
									 rx=rx, ry=ry,
									 stroke="black",
									 fill=all_chr_posi[chr_tmp]['chr'][3]))

	# 绘制chr name
	for chr_tmp in all_chr_label:
		dwg.add(dwg.text(all_chr_label[chr_tmp]["title"],
						 insert=(all_chr_label[chr_tmp]["x"], all_chr_label[chr_tmp]["y"]),
						 font_size="24",
						 fill="black",
						 stroke="none",
						 text_anchor="middle"))

	# 绘制marks
	for mark_tick in all_marks_ticks:
		dwg.add(dwg.line((mark_tick["x1"], mark_tick["y1"]), (mark_tick["x2"], mark_tick["y2"]), stroke=mark_tick["color"]))

	# 绘制label
	for mark_label in all_marks_labels:
		dwg.add(dwg.text(mark_label["title"], insert=(mark_label["x"], mark_label["y"] + 9),
						 font_size="18",
						 fill=mark_label["color"],
						 stroke=mark_label["color"]))

if __name__ == '__main__':
	conf = parse_config('./conf/draw.ini')
	log = gen_logger('./log/draw.log')

	dwg = svgwrite.Drawing(conf['global_settings']['output_filename'], profile='tiny',
						   size=(int(conf['image_settings']['image_width']),
								 int(conf['image_settings']['image_height']))
						   )
	dwg.add(dwg.text('Cotton Genome', insert=(1250, 40),
					 font_size="30",
					 fill="black",
					 stroke="none",
					 text_anchor="middle"))

	draw_pic(dwg, conf)

	dwg.save()

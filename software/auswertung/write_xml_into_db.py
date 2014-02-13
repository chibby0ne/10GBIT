#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright
# (C) 2011 creonic UG (haftungsbeschränkt)
# (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
#
# This file is part of the creonic simulation environment (CSE)
# for communication systems.
#
# \file
# \brief	XML file backend and result data handling
# \authors Timo Lehnigk-Emden, Matthias May, Christian Brehm
# \date	2011/12/20
#


from sim_db import *
from optparse import OptionParser 


if __name__ == "__main__":
	# parse command line arguments
	parser = OptionParser("write_xml_into_db.py [-d DATABASE] [-u USERNAME] XML_FILE_NAME ...") 
	parser.add_option("-d", dest="dbfile", default="")
	parser.add_option("-u", dest="username", default="")
	(opt, args) = parser.parse_args() 
	if len(args) < 1:
	    parser.error("XML file name missed.")
	if opt.username == '' or opt.dbfile == '':
		try:
			configs = open ('config.txt', 'r')
			f = configs.readlines()
		except:
			configs = open (os.path.dirname(__file__)+'/config.default', 'r')
			f = configs.readlines()
		if opt.dbfile == '':
			opt.dbfile = f[1][:len(f[1])-1]
		if opt.username == '':
			opt.username = f[3][:len(f[3])-1]

	tt = Sim_DB(opt.dbfile, opt.username)
	for a in args:
		tt.Process_New_File(str(a))
	tt.Show_DB()
	tt.Save_And_Close()


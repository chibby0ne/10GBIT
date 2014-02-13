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
# \authors Timo Lehnigk-Emden, Matthias May
# \date	2011/05/11
#


from sim_db import *
from optparse import OptionParser 


if __name__ == "__main__":
	# parse command line arguments
	parser = OptionParser("write_xml_into_db.py DATABASE") 
	(opt, args) = parser.parse_args() 
	if len(args)<1:
		parser.error("Database file name missed.")

	tt = Sim_DB(args[0], '')
	tt.Create_Tables()
	tt.Close()


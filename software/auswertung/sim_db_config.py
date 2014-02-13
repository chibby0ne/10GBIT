#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright
# (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
#
## \file
## \brief	Configuration for sim_db object.
#           Copy this file to sim_db_config.py and add the configuration for your chains here.<br><br>
## \author  Matthias May
#  \date	2011/10/27--2012/07/04
#


## Custom configuration class for Sim_DB object.
#
#  This file provides examples for configuring your Sim_DB object.
#  Copy this file to sim-db_config.py and add the configuration for your chains here.
class Sim_DB_Config(object):

	## Create this object and store the assigned Sim-DB object
	def __init__(self, sim_db_object):
		## Assigned Sim-DB object
		self.sdb = sim_db_object


	## Configure the assigned Sim-DB object
	def Config(self):
		sdb = self.sdb

		# Assign standards to chains
		#sdb.standards['mychain'] = ''

		## Names (printed in figure title) assigned to chains
		#sdb.chain_names['mychain'] = ''

		# Don't store these parameters in the configuration table:
		#sdb.param_blacklist['mychain'] = ["num_hi","il"]

		# List of chains with half iterations (instead of full iterations)
		#sdb.half_iter.append(
		#	'mychain')

		# Name of error rates module (default: Statistics_Error_Rates)
		#sdb.error_module['mychain'] = 'modulename'

		## Name of module and port providing the SNR; description for figure
		#
		# Default: module='global', port='eb_n0', description='E_b/N_0 / dB'
		# Format:
		#       self.snr[CHAIN] = [MODULE, PORT, DESCRIPTION]
		#sdb.snr['mychain'] = ['global', 'es_n0', 'E_s/N_0 / dB']

		# Name of iterating module and the port for number of iterations
		#sdb.iter_module_port['mychain'] = ['modulename', 'portname']

		# Customn result functions
		#                    name               function               module         port      y axis label  log scale (on/off)
		#sdb.result_function['sphere branch metrics'] = ('Get_Results_Over_Iterations', 'Mimo_Detection_Sphere',  'mean_num_bm_calcs', 'average number of branch metric calculations',       0)
		#sdb.result_function['sphere leaf nodes'] = ('Get_Results_Over_Iterations', 'Mimo_Detection_Sphere',  'mean_num_leaves', 'average number of leaf nodes',       0)
		#sdb.result_function['sphere nodes'] = ('Get_Results_Over_Iterations', 'Mimo_Detection_Sphere',  'mean_num_nodes', 'average number of visited nodes',       0)
		#sdb.result_function['sphere branch metrics hist'] = ('Get_Results_Over_Iterations', 'Mimo_Detection_Sphere',  'hist_num_bm_calcs', 'histogram of branch metric calculations',       0)
		#sdb.result_function['sphere leaf nodes hist'] = ('Get_Results_Over_Iterations', 'Mimo_Detection_Sphere',  'hist_num_leaves', 'histogram of leaf nodes',       0)
		#sdb.result_function['sphere nodes hist'] = ('Get_Results_Over_Iterations', 'Mimo_Detection_Sphere',  'hist_num_nodes', 'histogram of visited nodes',       0)
		#sdb.result_function['sphere error rates'] = ('Get_Results_Over_Iterations', 'error_rates_detection',  'error_rate_blocks', 'FER',       1)


	# Add, calculate, erase parameters (for title and legend of the plot) ....
	#
	# Set the format of a parameter by calling Append_Par().
	# The parameters will occur in the title/legend of the figure in the same order as Append_Par() is called in this function.
	# You can also add new parameters as, for instance, calculate the code rate.
	# The function Erase_Par() avoids printing of the given parameter.
	#
	# You can access the parameter values by calling Par() for arrays of string values,
	# Par_Int() or Par_Float() for arrays of integer or float values, repectively.
	# Float values can be converted to formatted strings by calling Format_Float().
	#
	# You can erase parameters under certain conditions by calling Mask() or Or(), respectively.
	# The functions Equal(), Equal_Single(), and Equal_Or_Default() compare the parameters to given values and return a truth mask.
	#
	# See sim_db_base_config.py for various examples.

	## Adjust parameters for a chain
	def Adjust_Param_mychain(self):
		sdb = self.sdb
		# example:
		#self.Append_Par('N', self.Par_Int('global.other_bits')+self.Par_Int('global.info_bits'))
		#self.Append_Par('K', self.Par('global.info_bits'))
		#self.Erase_Par('Channel_AWGN.start_seed')


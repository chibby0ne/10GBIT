#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright
# (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
#
## \file
## \authors Matthias May
## \brief	Basic configuration for sim_db object.
#           Add only the configuration for chains included in emssim.
#           For all other chains, add the configurations in sim_db_config.py.
#  \date    2011/10/27--2012/07/04
#


from numpy import *


## Basic configuration class for sim_db object.
#
#  Add only the configuration for chains included in emssim.
#  For all other chains, add the configurations in sim_db_config.py.
class Sim_DB_Base_Config(object):

	## Create this object and set configuration variables
	def __init__(self):

		## Standards assigned to chains
		self.standards = {}
		self.standards['emssim_lte_turbo'] = 'LTE'
		self.standards['emssim_hspa_turbo'] = 'HSPA'
		self.standards['emssim_umts_turbo'] = 'UMTS'

		## Names (printed in figure title) assigned to chains
		self.chain_names = {}
		self.chain_names['emssim_lte_turbo'] = 'LTE turbo decoder'
		self.chain_names['emssim_hspa_turbo'] = 'HSPA turbo decoder'
		self.chain_names['emssim_umts_turbo'] = 'UMTS turbo decoder'
		self.chain_names['rs_hard_kette'] = 'Reed Solomon Hard Decoder'

		## Parameter blacklist: Don't store these parameters in the configuration table.
		#
		# Note that module paramters taken from global parameters are automatically left out.
		self.param_blacklist = {}
		self.param_blacklist['emssim_lte_turbo'] = [
			"global.rec_poly","global.parity_poly","global.termination1","global.termination2","global.num_hi","global.il","global.bw_llrapp_out","global.noise_variance",
			"Statistics_Error_Rates.print_status_permanent","Statistics_Error_Rates.out_port_inner_dim_name","Statistics_Error_Rates.out_port_inner_dim_addr_offset",
			"Statistics_Error_Rates.max_num_diff_blocks","Statistics_Error_Rates.max_num_total_blocks",
			"Depunct_LTE.num_bits_after_depuncturing","Punct_LTE.num_bits_after_depuncturing"]
		self.param_blacklist['emssim_hspa_turbo'] = self.param_blacklist['emssim_lte_turbo']
		self.param_blacklist['emssim_umts_turbo'] = self.param_blacklist['emssim_lte_turbo']
		#self.param_blacklist['mimo_convolutional_chain_release'] = ["global.noise_variance", "global.es_n0"]
		self.param_blacklist['WPAN'] = ["global.noise_variance","global.es_n0"]


		## List of chains with half iterations (instead of full iterations)
		self.half_iter = [
			'emssim_lte_turbo',
			'emssim_hspa_turbo',
			'emssim_umts_turbo']

		## Name of error rates module (default: Statistics_Error_Rates)
		self.error_module = {}
		#self.error_module['mimo_convolutional_chain_release'] = 'error_rates_decoding'
		self.error_module['WPAN'] = 'error_rates_decoding'

		## Name of module and port providing the SNR; description for figure
		#
		# Default: module='global', port='eb_n0', description='E_b/N_0 / dB'
		# Format:
		#       self.snr[CHAIN] = [MODULE, PORT, DESCRIPTION]
		self.snr = {}
		#self.snr['mimo_convolutional_chain_release'] = ['global', 'eb_n0', 'E_b/N_0 / dB']
		self.snr['WPAN'] = ['global', 'eb_n0', 'E_b/N_0 / dB']

		## Name of iterating module and the port for number of iterations
		#
		# Format:
		#                      chain                 module     port
		self.iter_module_port = {}
		self.iter_module_port['emssim_lte_turbo'] = ['TCDec', 'hi_needed']
		self.iter_module_port['emssim_hspa_turbo'] = ['TCDec', 'hi_needed']
		self.iter_module_port['emssim_umts_turbo'] = ['TCDec', 'hi_needed']
		self.iter_module_port['WPAN'] = ['Decoder_LDPC_IEEE_802_15_3c', 'mean_iterations']


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


	## Adjust parameters for chain emssim_lte_turbo
	def Adjust_Param_emssim_lte_turbo(self):
		self.Append_Par('N', self.Par_Int('global.other_bits')+self.Par_Int('global.info_bits'))
		self.Append_Par('K', self.Par('global.info_bits'))
		try:
			self.Append_Par('R', self.Format_Float(self.Par_Float('global.info_bits')/(self.Par_Float('global.other_bits')+self.Par_Float('global.info_bits')), '%1.2f'))
		except:
			self.Append_Par('R', self.False_Array())
		self.Append_Par('Floating-point', self.Par('global.floating_point'))
		nofloat = self.Not(self.Par('global.floating_point'))
		self.Append_Par('Q_{\lambda}', self.Mask(self.Par('global.bw_in'), nofloat))
		self.Append_Par('Prec', self.Mask(self.Par('global.prec_in'), nofloat))
		self.Append_Par('Q_{\Lambda^e}', self.Mask(self.Par_Int('global.bw_in')+self.Par_Int('global.bw_llrextr_plus'), nofloat))
		self.Append_Par('Q_{\gamma}', self.Mask(self.Mask(self.Par('global.bw_gamma_sat'), nofloat), self.Not(self.Equal_Or_Default(self.Par('global.bw_gamma_sat'), '0'))))
		self.Append_Par('', self.Replace(self.Par('global.decoding_algo'), {
			"MAX_LOG_MAP_ESF":"Max-Log-MAP",
			"LOG_MAP":"Log-MAP",
			"MAX_LOG_MAP_ESF_BERROU_MM":"Max-Log-MAP, x+z saturation with one additional bit",
			"MAX_LOG_MAP_ESF_BERROU_TI":"Max-Log-MAP-ESF, x+z saturation",
			"MAX_LOG_MAP_ESF_BERROU_ASIP_OLD":"Max-Log-MAP-ESF, old ASIP, x+z saturation"}))
		self.Append_Par('ESF', self.Mask(self.Par('global.esf'), self.Not(self.Equal_Single(self.Par('global.decoding_algo'), 'LOG_MAP'))))
		print self.Par('global.decoding_algo')
		print self.Equal_Single(self.Par('global.decoding_algo'), 'LOG_MAP')
		print self.Not(self.Equal_Single(self.Par('global.decoding_algo'), 'LOG_MAP'))
		iswindowing = self.Not(self.Or(self.Equal(self.Par('global.info_bits'), self.Par('global.win_len')), self.Equal_Or_Default(self.Par('global.win_len'), '0')))
		self.Append_Par('', self.Mask(self.Par('global.map_architecture'), iswindowing))
		self.Append_Par('L_W', self.Mask(self.Par('global.win_len'), iswindowing))
		self.Append_Par('L_{ACQ}', self.Mask(self.Par('global.win_acq_len'), iswindowing))
		self.Append_Par('L_{sub-block ACQ}', self.Mask(self.Par('global.border_acq_len'), iswindowing))
		self.Append_Par('NII', self.Mask(self.Par('global.use_nii'), iswindowing))
		self.Append_Par('P_S', self.Mask(self.Par('global.parallel'), self.Mask(
			self.Not(self.Equal_Or_Default(self.Par('global.parallel'), '1')), self.Equal_Single(self.Par('global.map_architecture'), 'PSMAP'))))
		self.Append_Par('', self.Mask(
			self.Replace(self.Par('global.last_beta_calc'),{'BC_ACQ':'last beta from ACQ', 'BC_REC':'Last beta from recursion', 'BC_AUTO':'FALSE'}),
			self.Equal_Single(self.Par('global.map_architecture'), 'PSMAP')))
		self.Append_Par('alternating ACQ', self.Mask(self.Par('global.alternating_win_acq'), iswindowing))
		self.Append_Par('recursion patterns', self.Mask(self.Par('global.rec_patterns'), self.Mask(
			self.Not(self.Equal_Or_Default(self.Par('global.rec_patterns'), '2')), self.Equal_Single(self.Par('global.map_architecture'), 'XMAP'))))
		self.Append_Par('\Delta_k', self.Mask(self.Par('global.k_axis_shift'), self.Equal_Single(self.Par('global.map_architecture'), 'XMAP')))
		self.Append_Par('hard output iteration control #_{HI}', self.Mask(self.Par('global.hard_ic_num_hi'), self.Not(self.Equal_Or_Default(self.Par('global.hard_ic_num_hi'), '0'))))
		self.Append_Par('', self.Replace(self.Par('global.mapping'), {
			'MAP_BPSK':'BPSK',
			'MAP_QPSK':'QPSK',
			'MAP_16_QAM':'16-QAM',
			'MAP_64_QAM':'64-QAM',
			'MAP_256_QAM':'256-QAM'}))

		self.Erase_Par('Source_Bits.start_seed')
		self.Erase_Par('Source_Bits.mode')
		self.Erase_Par('Channel_AWGN.start_seed')
		self.Erase_Par('Demapper.channel_reliability')
		self.Erase_Par('global.crc_poly')
		self.Erase_Par('global.sm_init_val')


	## Adjust parameters for chain emssim_hspa_turbo
	def Adjust_Param_emssim_hspa_turbo(self):
		self.Adjust_Param_emssim_lte_turbo()


	## Adjust parameters for chain emssim_umts_turbo
	def Adjust_Param_emssim_umts_turbo(self):
		self.Adjust_Param_emssim_lte_turbo()
		
	## Adjust parameters for chain emssim_hspa_turbo
	def Adjust_Param_mimo_scenario_convolutional(self):
		self.Erase_Par('global.es_n0')


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
## \file
## \brief	XML file backend and result data handling
## \author	Timo Lehnigk-Emden, Matthias May
#  \date	2011/05/11--2012/07/04
#


from xml.dom import minidom
import logging
from logging import exception
import sqlite3
from time import *
import re  # regular expressions
from numpy import *
from sim_db_base_config import *
## Flag for avalability of sim_db_config.py
CONFIG = 1
try:
	from sim_db_config import *
except ImportError:
	CONFIG = 0


## Simulation database read and write access class.
#
# This class provides functions for
# \li creating a database for results
# \li reading XML result files and inserting the results into the database
# \li searching results in the database and catch them
# \li writing a Matlab file which plots the results
class Sim_DB(Sim_DB_Base_Config):

	## Create this object and open the database
	def __init__(self,  filename = ":memory:", username = "anonymous"):
		Sim_DB_Base_Config.__init__(self)

		# Open database
		## SQL database object
		self.db = sqlite3.connect(filename)
		## database username
		self.username = username

		## List of result functions
		self.result_function = {}
		self.result_function['FER over iterations'] = ('Get_Results_Over_Iterations', 'Error_Module()', 'error_rate_blocks', 'FER', 1)
		self.result_function['BER over iterations'] = ('Get_Results_Over_Iterations', 'Error_Module()', 'error_rate_bits', 'BER', 1)
		self.result_function['FER (single)'] = ('Get_Results_Single', 'Error_Module()', 'error_rate_blocks', 'FER', 1)
		self.result_function['BER (single)'] = ('Get_Results_Single', 'Error_Module()', 'error_rate_bits', 'BER', 1)
		self.result_function['Mean number of iterations'] = ('Get_Results_Single', 'Iter_Module()', 'Iter_Port()', "Mean number of iterations", 0)
		self.result_function['Mean number of iterations over max. iter.'] = ('Get_Results_Over_Iterations', 'Iter_Module()', 'Iter_Port()', "Mean number of iterations", 0)

		# Call custom configuration function
		if CONFIG:
			## Custom configuration object
			self.config = Sim_DB_Config(self)
			self.config.Config()
	
		# Set default values
		## Name of simulation chain
		self.chain = ''
		self.Reset_Results()



	# MAIN FUNCTIONS. Intended for calling in programs using this class. :::::::::::::::::::::::::::

	## Get simulations results
	def Get_Results(self, chain, result_type, names, values, hi_sel):
		## Current simulation chain
		self.chain = str(chain)
		print chain, result_type, self.result_function[str(result_type)]

		if str(result_type) in self.result_function:
			(function, module, port, ylabel, logscale) = self.result_function[str(result_type)]
			print function, module, port, ylabel, logscale

			if module=='Error_Module()':
				module = self.Error_Module()
			elif module=='Iter_Module()':
				module = self.Iter_Module()
			if port=='Error_Port()':
				port = self.Error_Port()
			elif port=='Iter_Port()':
				port = self.Iter_Port()
			if ylabel=='Mean number of iterations' and self.chain in self.half_iter:
				ylabel = 'Mean number of half-iterations'

			if getattr(self, function, 0):
				return getattr(self, function)(names, values, hi_sel, module, port, ylabel, logscale)

		print "Function not implemented."


	## Create tables
	def Create_Tables(self):

		# Create simulation table
		table_struct = "id INTEGER PRIMARY KEY ASC, \
		file_name TEXT, \
		file_content TEXT, \
		insert_date TEXT, \
		simulation_date TEXT, \
		user TEXT, \
		chain TEXT, \
		standard TEXT, \
		free_comment TEXT "
		self.db.execute("CREATE TABLE simulation ("+table_struct+")")
		
		# create configuration table
		table_struct = "id INTEGER PRIMARY KEY ASC, \
		module TEXT, \
		name TEXT, \
		value_float REAL, \
		value_string TEXT, \
		simulation_id INTEGER, \
		configuration_number INTEGER "
		self.db.execute("CREATE TABLE configuration ("+table_struct+")")

		# create result table
		table_struct = "id INTEGER PRIMARY KEY ASC, \
		module_name TEXT, \
		port_name TEXT, \
		inner_dim1_name TEXT, \
		inner_dim1_addr INTEGER, \
		inner_dim2_name TEXT, \
		inner_dim2_addr INTEGER, \
		value REAL, \
		simulation_id INTEGER, \
		configuration_number INTEGER "
		self.db.execute("CREATE TABLE result ("+table_struct+")")

		# create indices for faster search
		self.db.execute("CREATE INDEX simula_stand_ind ON simulation (standard);")
		self.db.execute("CREATE INDEX simula_chain_ind ON simulation (chain);")
		self.db.execute("CREATE INDEX config_simul_ind ON configuration (simulation_id);")
		self.db.execute("CREATE INDEX config_param_ind ON configuration (simulation_id,configuration_number);")
		self.db.execute("CREATE INDEX result_simul_ind ON result (simulation_id);")
		self.db.execute("CREATE INDEX result_param_ind ON result (simulation_id,configuration_number);")


	## Close connection to the database
	def Close(self):
		self.db.close()


	## Save tables and close connection to the database
	def Save_And_Close(self):
		self.db.commit()
		self.db.close()


	## Print all or the given row of the database
	def Show_DB(self, row_index = -1):

		if row_index == -1:
			print "simulation:"
			f = self.db.execute('SELECT id, file_name, "FILECONTENT NOT PRINTED", insert_date, simulation_date, user, chain, standard, free_comment FROM simulation')
			for row in f:
				print row

			#print "configuration:"
			#f = self.db.execute('SELECT * FROM configuration WHERE name<>"il"')
			#for row in f:
			#	print row
			#print "(interleaver not printed)"

			#print "result:"
			#f = self.db.execute('SELECT * FROM result')
			#for row in f:
			#	print row

		else:
			f = self.db.execute('SELECT * FROM simulation WHERE (rowid == ?)',  [row_index])
			for row in f:
				print row


	## Get all chains from those are results stored in the database
	def Get_Chains(self):
		chains = []
		f = self.db.execute('SELECT DISTINCT chain FROM simulation')
		for row in f:
			print str(row[0])
			chains.append(str(row[0]))
		return chains
		

	## Get all parameters used in the given chain
	def Get_Parameters(self, chain):
		params = []
		f = self.db.execute('SELECT DISTINCT module,name FROM configuration as c WHERE exists(select id from simulation where id=c.simulation_id and chain="%s")' % chain)
		for row in f:
			print str(row[0])+"."+str(row[1])
			params.append(str(row[0])+"."+str(row[1]))
		return params


	## Open an XML file and parse the data into the database
	def Process_New_File(self,  xml_filename):
		try:
			xmldoc = minidom.parse(xml_filename)
		except exception:
			print("File \"" + xml_filename + "\" cannot be opened!")

		print " NEW FILE: ",xml_filename, ' Try to read in!'
		self.Read_Save_Port_Data(xmldoc, xml_filename)


	## Delete a simulation in the database
	#
	# @param sim_id       identifier of the simulation to delete
	def Delete_Simulation(self, sim_id):
		print "Delete simulation with id=%d..." % int(sim_id)
		f = self.db.execute('DELETE FROM result where simulation_id = %d' % int(sim_id))
		f = self.db.execute('DELETE FROM configuration where simulation_id = %d' % int(sim_id))
		f = self.db.execute('DELETE FROM simulation where id = %d' % int(sim_id))


	## Write Matlab script for drawing figure
	def Write_Matlab_Script(self, filename):
		print "Writing Matlab script %s" % filename
		figfilename = re.sub('[{}\:\=,]','',self.title);
		figfilename = re.sub('[\ \/]','_',figfilename);
		fobj = open("%s" % filename, "w")
		fobj.write('clear all;\n')
		for i, sim in enumerate(self.xval):
			fobj.write('xval{%d} = [' % (i+1))
			for val in sim:
				fobj.write('%f ' % val)
			fobj.write('];\n')
		for i, sim in enumerate(self.yval):
			fobj.write('yval{%d} = [' % (i+1))
			for val in sim:
				fobj.write('%e ' % val)
			fobj.write('];\n')
		fobj.write('legendval = {')
		for val in self.legend:
			fobj.write('\'%s\' ' % val)
		fobj.write('};\n')
		fobj.write('graphtitle = \'%s\';\n' % self.title)
		fobj.write('figfilename = \'%s\';\n' % figfilename)
		fobj.write('xtitle = \'%s\';\n' % self.xtitle)
		fobj.write('ytitle = \'%s\';\n' % self.ytitle)
		fobj.write("%Farbvektor initialisieren, entsprechende Kurven haben gleiche Farben;\n\
		farb=zeros(1,3);\n\
		%Moegliche Marker\n\
		mark=['+','o','*','s','d','h'];\n\
		%Graph erzeugen\n\
		figf=figure('Visible','off');\n\
		for k=1:length(legendval)\n\
		    %Graph zeichnen, Farbe zufaellig, marker aus Vektor mark\n\
		    plot(xval{k}(:), yval{k}(:), 'linestyle','--', 'linewidth',1.2, 'color',rand(1,3), 'marker',mark(mod(k,6)+1), 'DisplayName',legendval{k});\n\
		    hold on;\n\
		end\n\
		%Graph formatieren\n\
		xlabel(xtitle);\n\
		ylabel(ytitle);\n\
		title(graphtitle);\n\
		axis tight;\n\
		grid on;\n")
		if self.logscale == 1:
			fobj.write("set(gca,'YScale','log');\n")
		fobj.write("set(gca,'YLimMode','auto');\n\
		set(gca,'Position',[0.05511 0.07292 0.9389 0.8847]);  %fits for create_pdf.m, paper size 24x18\n\
		box on;\n\
		hold off;\n\
		leg = legend('toggle');\n\
		legend('boxon');\n\
		set(leg,'Location','southwest','FontSize',10);\n\
		% Grafiken anzeigen\n\
		set(figf,'Visible','on');\n\
		% Grafik speichern\n\
		%saveas(figf,[figfilename, '.fig'],'fig');\n\
		%h = gcf();\n\
		%set(h,'PaperSize',[24 18]);\n\
		%set(h,'PaperUnits','centimeters');\n\
		%set(h,'PaperOrientation','portrait');\n\
		%set(h,'PaperPosition',[0 0 24 18]);\n\
		%print('-dpdf',[figfilename, '.pdf']);\n")
		fobj.close()



	# HELPER FUNCTIONS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	## Get name of error rates module
	def Error_Module(self):
		if self.chain in self.error_module:
			return self.error_module[self.chain]
		return 'Statistics_Error_Rates'


	## Get name of the iterating module in the chain
	def Iter_Module(self):
		if self.chain in self.iter_module_port:
			return self.iter_module_port[self.chain][0]
		exit("No iterating module for the chain "+self.chain+" specified in the configuration")


	## Get the iteration number port of the iterating module in the chain
	def Iter_Port(self):
		if self.chain in self.iter_module_port:
			return self.iter_module_port[self.chain][1]
		exit("No iterating module for the chain "+self.chain+" specified in the configuration")


	## Get SNR Module
	def Snr_Module(self):
		if self.chain in self.snr:
			return self.snr[self.chain][0]
		else:
			return 'global'


	## Get SNR Port
	def Snr_Port(self):
		if self.chain in self.snr:
			return self.snr[self.chain][1]
		else:
			return 'eb_n0'


	## Get SNR Description
	def Snr_Desc(self):
		if self.chain in self.snr:
			return self.snr[self.chain][2]
		else:
			return 'E_b/N_0 / dB'


	## Get SNR Module.Port
	def Snr_Module_Port(self):
		return self.Snr_Module()+"."+self.Snr_Port()


	## Get the type of a value (as column name in table configuration)
	def Get_Value_Col(self, value):
		try:
			value_typed = float(value)
			try:
				if 'x' in value or 'X' in value:
					value_col = 'value_string'
				else:
					value_col = 'value_float'
			except:
				value_col = 'value_float'
		except:
			value_col = 'value_string'
		return value_col


	## Get the type of the first value in the list (as column name in table configuration)
	def Get_Value_List_Col(self, value):
		val_list = eval(value)
		try:
			val0 = val_list[0]
		except:
			val0 = val_list  # only one entry: val is not a list
		return self.Get_Value_Col(val0)


	## Get a value returned as float or string			   
	def Get_Value_Typed(self, value):
		try:
			value_typed = float(value)
		except:
			value_typed = str(value).upper()
		return value_typed


	## Write a configuration entry
	def Write_Config_Entry(self,  configuration_number, module, name, value, auto_param_blacklist):
		if not ((self.chain in self.param_blacklist and (str(module)+"."+str(name)) in self.param_blacklist[self.chain])
			or ((str(module)+"."+str(name)) in auto_param_blacklist)):
			insert_str = 'INSERT INTO configuration (module, name, '+self.Get_Value_Col(value)+', simulation_id, configuration_number) VALUES (?, ?, ?, ?, ?)'
			data = []
			data.append(str(module))
			data.append(str(name))
			data.append(self.Get_Value_Typed(value))
			data.append(self.simulation_id)
			data.append(int(configuration_number))
			self.db.execute(insert_str, data)


	## Read the data out of the XML tree and store them into a database table
	def Read_Save_Port_Data(self,  xml, filename):

		# Read simulation data
		data = []
		data.append(filename)
		# get complete content of result file							
		fobj = open(filename, "r")
		data.append(fobj.read())
		fobj.close()
		insert_date = localtime()
		data.append(strftime("%c", insert_date))
		simulation_date = ""
		if len(xml.getElementsByTagName("simulation_date"))>0:
			simulation_date = str(xml.getElementsByTagName("simulation_date")[0].firstChild.nodeValue)
		else:
			print "simulation_date NOT FOUND"
		data.append(simulation_date)
		data.append(self.username)
		exec_name = ""
		if len(xml.getElementsByTagName("executable_name"))>0:
			exec_name = re.search('[^/]+\Z', str(xml.getElementsByTagName("executable_name")[0].firstChild.nodeValue)).group(0)
			exec_name = re.sub('\.[a-zA-Z0-9]+','',exec_name)
		else:
			print "exec_name NOT FOUND"
		data.append(exec_name)
		if exec_name in self.standards:
			data.append(self.standards[exec_name])
		else:
			data.append('none')
		if len(xml.getElementsByTagName("free_comment"))>0:
			free_comment = str(xml.getElementsByTagName("free_comment")[0].firstChild.nodeValue)
		else:
			free_comment = ""
		data.append(free_comment)
		insert_str = "INSERT INTO simulation (file_name, file_content, insert_date, simulation_date, user, chain, standard, free_comment) VALUES (?, ?, ? , ?, ?, ?, ?, ?)"
		self.db.execute(insert_str,  data)
		self.chain = exec_name

		# Read simulation id (which the database has assigned to the current simulation)
		query = "SELECT id FROM simulation WHERE chain = '%s' AND simulation_date = '%s'" % (exec_name, simulation_date)
		f = self.db.execute(query)
		## Current simulation id
		self.simulation_id = 0
		for row in f:
			self.simulation_id = row[0]

		# Collect all global parameters used in modules (do not store them later)
		auto_param_blacklist = []
		init_config = xml.getElementsByTagName("initial_configuration")[0]
		for module in init_config.getElementsByTagName("module"):
			module_name = str(module.getElementsByTagName("instance_name")[0].firstChild.nodeValue)
			for parameter in module.childNodes:
				name = parameter.localName
				if name and name != "instance_name":
					if len(parameter.getElementsByTagName("global_variable")) > 0:
						auto_param_blacklist.append(str(module_name)+"."+str(name))

		# Gor for all configuration sections
		config = xml.getElementsByTagName("configuration")
		for i, item in enumerate(config):

			# Extract all global parameters
			glob = item.getElementsByTagName("global")
			for j, jtem in enumerate(glob):
				parameters = jtem.getElementsByTagName("variable")
				for k, ktem in enumerate(parameters):
					self.Write_Config_Entry(i, 'global', ktem.getAttribute("name"), ktem.firstChild.nodeValue, auto_param_blacklist)

			# Extract all module parameters
			module = item.getElementsByTagName("module")
			for j, jtem in enumerate(module):
				module_name = str(jtem.getElementsByTagName("instance_name")[0].firstChild.nodeValue)
				parameters = jtem.childNodes
				for k, ktem in enumerate(parameters):
					name = ktem.localName
					if name and name != "instance_name":
						self.Write_Config_Entry(i, module_name, name, ktem.firstChild.nodeValue, auto_param_blacklist)

		# Go for all result sections
		result = xml.getElementsByTagName("result")
		for i,  item in enumerate(result):
			#print i

			# Extract all modules
			module = item.getElementsByTagName("module")
			for j, jtem in enumerate(module):
				module_name = jtem.getAttribute("name")
				#print "  ", j, module_name

				# Extract all status_out port of each module
				port = jtem.getElementsByTagName("status_out")
				for k,  ktem in enumerate(port):
					port_name = ktem.getAttribute("name")
					subval = ktem.getElementsByTagName("value")
					#print "     ", k, port_name, subval

					# Extract all values in the port
					for m,  mtem in enumerate(subval):
						raw_val = float(mtem.firstChild.nodeValue)
						#print "        ", m, raw_val

						# Get first attribute (inner dimension name and value)
						attr = mtem.attributes.item(0)
						#print "           attr=", attr

						if(attr!=None):
							dim_name = attr.name
							dim_value = float(attr.value)
							insert_str = "INSERT INTO result (module_name,port_name,inner_dim1_name,inner_dim1_addr,value,configuration_number,simulation_id) VALUES (?,?,?,?,?,?,?)"
							data = []
							data.append(module_name)
							data.append(port_name)
							data.append(dim_name)
							data.append(dim_value)
							data.append(raw_val)
							data.append(i)
							data.append(self.simulation_id)
						else:
							insert_str = "INSERT INTO result (module_name,port_name,value,configuration_number,simulation_id) VALUES (?,?,?,?,?)"
							data = []
							data.append(module_name)
							data.append(port_name)
							data.append(raw_val)
							data.append(i)
							data.append(self.simulation_id)

						# insert into result table
						#print 'RESULT:', insert_str, data
						self.db.execute(insert_str,  data)


	## Returns all indices of all items in the list l whose value is x
	def Get_All(self, x, l):
		ind = []
		for i, li in enumerate(l):
			if li==x:
				ind.append(i)
		return ind


	## Returns cut-set of two lists
	def Get_Cut(self, a, b):
		c = []
		for e in a:
			if e in b:
				c.append(e)
		return c


	## Returns 1 if all elements in the list are indentical, otherwise 0
	def Check_All_Identical(self, l):
		if len(l)>0:
			for a in l[1:]:
				if a != l[0]:
					return 0
		return 1


	## Get Value as a string; and without .0, if integer
	def Get_Value_Human(self, value_float, value_string):
		if value_float != None:
			# trim x.0 to x
			if (int(value_float) - value_float) == 0:
				val_human = str(int(value_float))
			else:
				val_human = str(value_float)
		else:
			val_human = str(value_string)
		return val_human


	## Convert a list to a string (for SQL queries like "WHERE x in (a,b,c)")
	def List_To_String(self, l):
		s = re.sub('\[', '(', str(l))
		s = re.sub('\]', ')', s)
		return s


	## Write Legend and Title
	def Write_Legend_Title(self, resultname, logscale):

		# Create legend (with iterations)
		if len(self.iter)>0:
			for h in self.iter:
				if h >= 0:
					if self.chain in self.half_iter:
						if h%2:
							h = "%1.1f" % (float(h)/2)
						else:
							h = "%d" % (h/2)
					self.legend.append("iteration %s" % h)
				else:
					self.legend.append("")
		else:
			self.legend.append("")

		# Adjust and insert parameters (chain specific)
		self.Adjust_Param()

		# Insert all remaining parameters
		for p in self.pardict:
			if not p in self.par_used:
				self.Append_Par(p, self.pardict[p])

		# Write axis titles and scale
		self.xtitle = self.Snr_Desc()
		self.ytitle = resultname
		self.logscale = logscale

		# Console output of found simulation results
		self.Print_Results()


	## Add, calculate, delete parameters (for title and legend of the plot)
	def Adjust_Param(self):
		if getattr(self, "Adjust_Param_"+self.chain, 0):
			getattr(self, "Adjust_Param_"+self.chain)()

		if CONFIG:
			if getattr(self.config, "Adjust_Param_"+self.chain, 0):
				getattr(self.config, "Adjust_Param_"+self.chain)()


	## Reset result lists
	def Reset_Results(self):
		## Figure title
		self.title = self.chain
		if self.chain in self.chain_names:
			self.title = self.chain_names[self.chain]
		## Figure legend
		self.legend = []
		## List of x values in figure
		self.xval = []
		## List of y values in figure
		self.yval = []
		## x axis title in figure
		self.xtitle = ''
		## y axis title in figure
		self.ytitle = ''
		## Flag for logarithmic scaling of y axis in figure
		self.logscale = 1
		## Dictionary for parameters and their values
		self.pardict = {}
		## List of already used parameters
		self.par_used = {}
		## List of chosen iterations
		self.iter = []
		## Number of iterations to select
		self.num_iter = 0
		## Number of added different simulations
		self.sims_added = 0
		## Current position of added data within one simulation
		self.data_index = [0]
		

	## Search configurations (simulation_id, configuration_id)
	def Search_Configurations(self, names, values):
		self.Reset_Results()
		query = 'SELECT DISTINCT simulation_id,configuration_number FROM configuration as c where exists (select id from simulation where id=c.simulation_id and chain="%s")' % self.chain
		# .. first all parameters where we search for default values...
		for n, pname in enumerate(names):
			mod = re.sub("\..*", "", pname)
			name = re.sub(".*\.", "", pname)
			val = str(values[n]).upper()
			if re.search("!", val):
				# search for default value (parameter not stored)
				query += ' and ((not exists (SELECT id FROM configuration WHERE simulation_id=c.simulation_id and configuration_number=c.configuration_number and name="%s"))' % name
				if re.search("\(!,", val):        
					val = re.sub("\(!,", "(", val)
					query += ' or (exists (SELECT id FROM configuration WHERE simulation_id=c.simulation_id and configuration_number=c.configuration_number and\
						module="%s" and name="%s" and %s in %s))' % (mod, name, self.Get_Value_List_Col(val), val)
				query += ')';
		# .. then the other parameters (as intersect query)...
		for n, pname in enumerate(names):
			mod = re.sub("\..*", "", pname)
			name = re.sub(".*\.", "", pname)
			val = str(values[n]).upper()
			if not re.search("!", val):
				if re.search("^<>", val):
					val = re.sub("<>", "", val)
					query += ' INTERSECT SELECT simulation_id,configuration_number FROM configuration WHERE module="%s" and name="%s" and %s <> \'%s\'' % (
						mod, name, self.Get_Value_Col(val), val)                                                                                                                # not equal
				elif re.search("^<", val):
					val = re.sub("<", "", val)
					query += ' INTERSECT SELECT simulation_id,configuration_number FROM configuration WHERE module="%s" and name="%s" and value_float < %s' % (mod, name, val)  # less than
				elif re.search("^>", val):
					val = re.sub(">", "", val)
					query += ' INTERSECT SELECT simulation_id,configuration_number FROM configuration WHERE module="%s" and name="%s" and value_float > %s' % (mod, name, val)  # greater than
				elif re.search("^\(.*\)", val):
					query += ' INTERSECT SELECT simulation_id,configuration_number FROM configuration WHERE module="%s" and name="%s" and %s in %s' % (
						mod, name, self.Get_Value_List_Col(val), val)                                                                                                           # list of values
				else:
					query += ' INTERSECT SELECT simulation_id,configuration_number FROM configuration WHERE module="%s" and name="%s" and %s=\'%s\'' % (
						mod, name, self.Get_Value_Col(val), val)                                                                                                                # just one value
		print "query:",query
		return self.db.execute(query)


	## Fetch module parameters
	def Fetch_Module_Param(self, config):
		parfetch = []
		query = 'SELECT module, name, value_float, value_string FROM configuration WHERE simulation_id=%d and configuration_number=%d' % (config[0], config[1])
		query += ' and not (module="%s" and name="%s")' % (self.Snr_Module(), self.Snr_Port())
		fp = self.db.execute(query)
		for row in fp:
			parfetch.append((str(row[0])+"."+str(row[1]), self.Get_Value_Human(row[2], row[3])))
		return parfetch


	## Fetch SNR value
	def Fetch_SNR(self, config):
		fx = self.db.execute('SELECT value_float FROM configuration WHERE simulation_id=%d and configuration_number=%d and module="%s" and name="%s"' % (
			config[0], config[1], self.Snr_Module(), self.Snr_Port()))
		snr = 0
		for row in fx:
			snr = float(row[0])
		return snr


	## Insert x and y values into lists
	def Insert_Result(self, yfetch, config):
		simid = -1;
	
		# Fetch parameters and SNR values
		parfetch = self.Fetch_Module_Param(config)
		xfetch = self.Fetch_SNR(config)

		# gucken, ob parameter ueberhaupt schonmal da war
		new_sim = 0
		for row in parfetch:
			if not row[0] in self.pardict:
				new_sim = 1
				print "parameter",row[0],"noch nicht vorhanden"
				# parameter ergaenzen, wert bei allen vorher auf default setzen
				vals = []
				for i in range(0,self.sims_added):
					vals.append("DEFAULT")
				self.pardict[row[0]] = vals

		if new_sim==0:
			# gucken, ob ein uebereinstimmender parametersatz vorhanden ist
			possible_sims = arange(0, self.sims_added)
			for row in parfetch:
				possible_sims = self.Get_Cut(possible_sims, self.Get_All(row[1], self.pardict[row[0]]))
			if len(possible_sims)==0:
				new_sim = 1
				print "parameter abweichend"
			else:
				simid = possible_sims[0]

		if new_sim==1:
			simid = self.sims_added
			self.sims_added += 1
			# add parameter value to dict
			for row in parfetch:
				#print "    append", row[0], row[1]
				self.pardict[row[0]].append(row[1])
			# write x and y values
			for r, row in enumerate(yfetch):
				self.iter.append(int(row[0]))
				self.yval.append([float(row[1])])
				self.xval.append([xfetch])
			if simid>0:
				self.data_index.append(self.data_index[simid-1]+len(yfetch))
			self.num_iter = len(yfetch)
		else:
			# write x and y values
			for r, row in enumerate(yfetch):
				ind = self.data_index[simid]+r
				toinsert = 1
				for i, val in enumerate(self.xval[ind]):
					if toinsert and xfetch <= val:
						self.yval[ind].insert(i, float(row[1]))
						self.xval[ind].insert(i, xfetch)
						toinsert = 0
				if toinsert:
					self.yval[ind].append(float(row[1]))
					self.xval[ind].append(xfetch)

		# set values of not available parameters (in current fetch) to DEFAULT
		for p in self.pardict.items():
			if len(p[1]) < self.sims_added:
				p[1].append("DEFAULT")


	## Return a string in the form /par/=/val/
	def Parstring(self, par, val):
		if str(val) == 'FALSE':
			return ''
		elif str(val) == 'TRUE':
			return par
		else:
			if par == '':
				return str(val)
			else:
				return par+"="+str(val)
		

	## Append a parameter to the title or the legend
	def Append_Par(self, par, val):
		if len(val) > 0:
			if self.Check_All_Identical(val)==1:
				# insert parameter into title
				if self.Parstring(par, val[0]) != '':
					self.title += ", "+self.Parstring(par, val[0])
			else:
				# insert parameter into legend
				for i, v in enumerate(val):
					if self.Parstring(par, v) != '':
						for h in range(self.num_iter):
							self.legend[i*self.num_iter+h] += str(", "+self.Parstring(par, v))


	## Mark the parameter as used (such that it is left out in the title/legend)
	def Erase_Par(self, par):
		self.par_used[par] = self.Val_Array(1)


	## Get the values of the given parameter as an integer array
	def Par_Int(self, par):
		if par in self.pardict:
			self.par_used[par] = ones( (len(self.pardict.itervalues().next())))
			return array(self.pardict[par], int)
		return self.Val_Array(0)


	## Get the values of the given parameter as an integer array
	def Par_Float(self, par):
		if par in self.pardict:
			self.par_used[par] = ones( (len(self.pardict.itervalues().next())))
			return array(self.pardict[par], float)
		return self.Val_Array(0.0)


	## Get the values of the given parameter
	def Par(self, par):
		if par in self.pardict:
			self.par_used[par] = ones( (len(self.pardict.itervalues().next())))
			return self.pardict[par]
		return self.False_Array()


	## Format the elements in the array val
	def Format_Float(self, val, form):
		try:
			for vi, v in enumerate(val):
				val[vi] = form % float(v)
			return val
		except:
			print 'ERROR: Format_Float: val =', val
			return self.False_Array()


	## Return a copy of the array val with the elements replaced according to the dictionary repl
	def Replace(self, val, repl):
		try:
			valnew = []
			for v in val:
				if str(v) in repl:
					valnew.append(repl[str(v)])
				else:
					valnew.append(v)
			return valnew
		except:
			print 'ERROR: Replace: val =', val
			return self.False_Array()


	## Returns a mask with 'TRUE' at the positions where the value in the array val1 is equal to the value in the array val2, else 'FALSE'
	def Equal(self, val1, val2):
		try:
			mask = []
			for vi, v1 in enumerate(val1):
				if v1 == val2[vi].upper():
					mask.append('TRUE')
				else:
					mask.append('FALSE')			
			return mask
		except:
			print 'ERROR: Equal: val1 =', val1, "  val2 =", val2
			return self.False_Array()


	## Returns a mask with 'TRUE' at the positions where the value in the list val is equal to comp, else 'FALSE'
	def Equal_Single(self, val, comp):
		try:
			mask = []
			for v in val:
				if v == comp.upper():
					mask.append('TRUE')
				else:
					mask.append('FALSE')			
			return mask
		except:
			print 'ERROR: Equal_Single: val =', val
			return self.False_Array()


	## Returns a mask with 'TRUE' at the positions where the value in val is equal to comp or 'DEFAULT', else 'FALSE'
	def Equal_Or_Default(self, val, comp):
		try:
			mask = []
			for v in val:
				if v == comp.upper() or v == 'DEFAULT':
					mask.append('TRUE')
				else:
					mask.append('FALSE')			
			return mask
		except:
			print 'ERROR: Equal_Or_Default: val =', val
			return self.False_Array()


	## Calculate OR of parameter values (if both are not 'FALSE', take value from val1)
	def Or(self, val1, val2):
		vala = []
		try:
			for vi, v1 in enumerate(val1):
				if v1.upper() == 'FALSE':
					vala.append(val2[vi])
				else:
					vala.append(v1)
			return vala
		except:
			print 'ERROR: Or: val1=', val1, "  val2=", val2
			return self.False_Array()


	## Invert parameter values ('TRUE' if mask='FALSE', else 'FALSE')
	def Not(self, val):
		newval = []
		for v in val:
			if v.upper() == 'FALSE':
				newval.append('TRUE')
			else:
				newval.append('FALSE')
		return newval


	## Mask parameter values ('FALSE' if mask='FALSE', else value)
	def Mask(self, val, mask):
		newpar = []
		try:
			for vi, v in enumerate(val):
				if mask[vi].upper() == 'FALSE':
					newpar.append('FALSE')
				else:
					newpar.append(v)
			return newpar
		except:
			return self.False_Array()


	## Returns an parameter value array with all entries = val
	def Val_Array(self, val):
		try:
			vala = []
			for i in arange(len(self.pardict.itervalues().next())):
				vala.append(val)
			return vala
		except:
			return []


	## Returns an parameter value array with all entries = 'FALSE'
	def False_Array(self):
		return self.Val_Array('FALSE')


	## Print found results on the console
	def Print_Results(self):
		print "xval:",self.xval
		print "yval:",self.yval
		print "title:",self.title
		print "legend:",self.legend


	## Get results over iterations
	#
	# @param names        list of parameter names
	# @param values       list of parameter values
	# @param iter_sel     list of iterations to select
	# @param module       name of the module which provides the results
	# @param port         name of the module's port which provides the results
	# @param resultname   result name, written at y-axis
	# @param logscale     1: logarithmic scaling for y-axis, 0: linear scaling for y-axis
	def Get_Results_Over_Iterations(self, names, values, iter_sel, module, port, resultname, logscale):
		print "\nSearching simulation results: %s over iterations:" % resultname, self.chain, names, values, iter_sel

		# fetch parameters and results for all found configurations; group if only SNR differs
		for c, config in enumerate(self.Search_Configurations(names, values)):
			print c, config			

			# read results for all iterations
			if len(iter_sel) > 0:
				fy = self.db.execute('SELECT inner_dim1_addr,value FROM result WHERE simulation_id=%d and configuration_number=%d and \
					module_name="%s" and port_name="%s" and inner_dim1_addr in %s' % (config[0], config[1], module, port, self.List_To_String(iter_sel)))
			else:
				fy = self.db.execute('SELECT inner_dim1_addr,value FROM result WHERE simulation_id=%d and configuration_number=%d and \
					module_name="%s" and port_name="%s"' % (config[0], config[1], module, port))
			yfetch = []
			for row in fy:
				yfetch.append(row)	

			self.Insert_Result(yfetch, config)

		# Write legend and title
		self.Write_Legend_Title(resultname, logscale)


	## Get results (single value per SNR point)
	#
	# @param names        list of parameter names
	# @param values       list of parameter values
	# @param iter_sel     IGNORED
	# @param module       name of the module which provides the results
	# @param port         name of the module's port which provides the results
	# @param resultname   result name, written at y-axis
	# @param logscale     1: logarithmic scaling for y-axis, 0: linear scaling for y-axis
	def Get_Results_Single(self, names, values, iter_sel, module, port, resultname, logscale):
		print "\nSearching simulation results: %s:" % resultname, self.chain, names, values

		# fetch parameters and results for all found configurations; group if only SNR differs
		for c, config in enumerate(self.Search_Configurations(names, values)):
			print c, config			

			# read result
			fy = self.db.execute('SELECT value FROM result WHERE simulation_id=%d and configuration_number=%d and \
				module_name="%s" and port_name="%s"' % (config[0], config[1], module, port))
			yfetch = []
			for row in fy:
				yfetch.append([-1, row[0]])  # put result at position; iter=-1 means no difference between iterations

			self.Insert_Result(yfetch, config)

		# Write legend and title
		self.Write_Legend_Title(resultname, logscale)


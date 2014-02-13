#!/usr/bin/python
# -*- coding: utf8 -*-
#
# Copyright
# (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
#
## \file
## \brief	Graphical user interface for getting simulation results out of the database.
## \author  Norbert Bass, Matthias May, Christian Brehm
# \date	    2011/12/20
#


import sys
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QFileDialog
from Ui import Ui_Dialog
from sim_db import * #import read functions
import xml.dom.minidom as dom
import xml.dom
import os

## Graphical user interface for getting simulation results out of the database
class Gui(QtGui.QDialog,Ui_Dialog):	
	
	## Create this object and login to the database
	def __init__(self, parent=None):
		QtGui.QDialog.__init__(self)
		self.setupUi(self)
		QtCore.QObject.connect(self.pushButton_reLogin, QtCore.SIGNAL("pressed()"), self.reLogin)
		QtCore.QObject.connect(self.comboBox_4, QtCore.SIGNAL("currentIndexChanged(QString)"), self.readParameter)
		QtCore.QObject.connect(self.pushButton_add, QtCore.SIGNAL("pressed()"), self.addIteration)
		QtCore.QObject.connect(self.pushButton_del_hi_sel, QtCore.SIGNAL("pressed()"), self.del_hi_sel)
		QtCore.QObject.connect(self.pushButton_search, QtCore.SIGNAL("pressed()"), self.search)
		QtCore.QObject.connect(self.pushButton_save, QtCore.SIGNAL("pressed()"), self.save)
		QtCore.QObject.connect(self.pushButton_saveAs, QtCore.SIGNAL("pressed()"), self.saveAs)
		QtCore.QObject.connect(self.pushButton_exit, QtCore.SIGNAL("pressed()"), self.exit)

		self.login()
		
		## Additional iterations to select
		self.hi_sel_help = []

	
	## Initial database login
	def login(self):
		# read configuration file
		try:
			configs = open ('config.txt', 'r')
			f = configs.readlines()
		except:
			print '\nFile config.txt cannot be opened. Reading default config file...'
			configs = open (os.path.dirname(__file__)+'/config.default', 'r')
			f = configs.readlines()
			self.writeConfigFile(f)
		
		database_fileName = f[1][:len(f[1])-1]
		self.lineEdit_db.clear()
		self.lineEdit_db.insert(database_fileName)
		usernameName = f[3][:len(f[3])-1]
		self.lineEdit_6.clear()
		self.lineEdit_6.insert(usernameName)
		configs.close()

		# try to open database
		try:
			fo = open (database_fileName, 'r')
			## Sim_Db object for access to the database
			self.tt = Sim_DB(database_fileName, usernameName)
		except:
			self.tt = 0
			
		if not (self.tt==0):
			for item in sort(self.tt.result_function.keys()):
				self.comboBox_5.addItem(item)
				self.comboBox_5.setCurrentIndex(int(f[9]));
		self.tableWidget.setColumnWidth(0,350)
		self.tableWidget.setColumnWidth(1,150)
		self.readChains()

	## Database login after changing the database
	def reLogin(self):
		if not (self.tt==0):
			self.tt.Close()
			del self.tt
		database_fileName = str(self.lineEdit_db.text())
		usernameName = str(self.lineEdit_6.text())
		try:
			fo = open (database_fileName, 'r')
			self.tt = Sim_DB(database_fileName, usernameName)
		except:
			self.tt = 0
		self.readChains()
	
	## Read alls chains for which the database contains simulation results
	def readChains(self):
		if self.tt==0:
			return
		self.comboBox_4.clear()
		chains = self.tt.Get_Chains()
		for c in chains:
			self.comboBox_4.addItem(str(c))
		configs = open ('config.txt', 'r')
		f = configs.readlines()
		self.comboBox_4.setCurrentIndex(int(f[7]));
	
	## Read all simulation parameters which are stored for the current chain
	def readParameter(self):
		if self.tt==0:
			return
		chain = self.comboBox_4.currentText()
		names = self.tt.Get_Parameters(chain)
		pc = ( len(names) )# parametercount
		self.tableWidget.setRowCount(pc)
		for i in range(pc):
			self.tableWidget.setItem(i, 0, QtGui.QTableWidgetItem(names[i]))
			self.tableWidget.setItem(i, 1, QtGui.QTableWidgetItem(''))
		configs = open ('config.txt', 'r')
		f = configs.readlines()
		configs.close()
		temp = f[5]
		x = len(temp)-1
		fileName = temp[:x]
		self.lineEdit_4.clear()
		self.lineEdit_4.insert(fileName)
	
	## Add an iteration to the iteration list for searching results
	def addIteration(self):
		self.hi_sel_help.append(self.iterationBox.value())
		self.listWidget_2.clear()
		self.hi_sel_help = list(set(self.hi_sel_help))
		self.listWidget_2.addItem(str(sort(self.hi_sel_help)))		

	## Clear iteration list for searching results
	def del_hi_sel(self):
		del self.hi_sel_help[:]
		self.listWidget_2.clear()

	## Search results in database (calls Sim_DB:GetResults())
	def search(self):
		if self.tt==0:
			return
		chain = self.comboBox_4.currentText()
		names = self.tt.Get_Parameters(chain)
		pc = len(names) #parametercount
		values = [None]*pc
		i = 0
		j = 0
		k = 0
		while i < pc:
			if self.tableWidget.item(i, 1) is None:
				del names[i-j]
				j += 1
			elif not str(self.tableWidget.item(i, 1).text()):
				del names[i-j]
				j += 1
			else:
				temp = self.tableWidget.item(i, 1).text()
				values[k] = str(temp)
				k += 1
			i += 1
		del values[k:]
		chain = self.comboBox_4.currentText()
		result_type = self.comboBox_5.currentText()
		hi_sel = []
		if self.checkBox.isChecked():
			hi_sel.append(0)
		if self.checkBox_1.isChecked():
			hi_sel.append(1)
		if self.checkBox_2.isChecked():
			hi_sel.append(2)
		if self.checkBox_3.isChecked():
			hi_sel.append(3)
		if self.checkBox_4.isChecked():
			hi_sel.append(4)
		if self.checkBox_5.isChecked():
			hi_sel.append(5)
		if self.checkBox_6.isChecked():
			hi_sel.append(6)
		if self.checkBox_7.isChecked():
			hi_sel.append(7)
		if self.checkBox_8.isChecked():
			hi_sel.append(8)
		if self.checkBox_9.isChecked():
			hi_sel.append(9)
		if self.checkBox_10.isChecked():
			hi_sel.append(10)
		if self.checkBox_11.isChecked():
			hi_sel.append(11)
		if self.checkBox_12.isChecked():
			hi_sel.append(12)
		if self.checkBox_13.isChecked():
			hi_sel.append(13)
		if self.checkBox_14.isChecked():
			hi_sel.append(14)
		if self.checkBox_15.isChecked():
			hi_sel.append(15)
		if self.checkBox_16.isChecked():
			hi_sel.append(16)
		if self.checkBox_17.isChecked():
			hi_sel.append(17)
		if self.checkBox_18.isChecked():
			hi_sel.append(18)
		if self.checkBox_19.isChecked():
			hi_sel.append(19)
		if self.checkBox_20.isChecked():
			hi_sel.append(20)
		if self.checkBox_21.isChecked():
			hi_sel.append(21)
		if self.checkBox_22.isChecked():
			hi_sel.append(22)
		if self.checkBox_23.isChecked():
			hi_sel.append(23)
		if self.checkBox_24.isChecked():
			hi_sel.append(24)
		if self.checkBox_25.isChecked():
			hi_sel.append(25)
		if self.checkBox_26.isChecked():
			hi_sel.append(26)
		if self.checkBox_27.isChecked():
			hi_sel.append(27)
		if self.checkBox_28.isChecked():
			hi_sel.append(28)
		if self.checkBox_29.isChecked():
			hi_sel.append(29)
		if self.checkBox_30.isChecked():
			hi_sel.append(30)
		if self.checkBox_31.isChecked():
			hi_sel.append(31)
		
		i = 0
		y = len(self.hi_sel_help)-1
		while i <= y:
			hi_sel.append(self.hi_sel_help[i])
			i += 1
		hi_sel = list(set(hi_sel))
		
		self.tt.Get_Results(chain, result_type, names, values, hi_sel)
		
		self.listWidget.clear()
		#self.listWidget.addItem("Title:")
		self.listWidget.addItem(self.tt.title)
		#self.listWidget.addItem("Legend:")
		lc = ( len(self.tt.legend) ) # legendcount
		i = 0
		while i < lc:
			self.listWidget.addItem(str(i+1) + ". " + (self.tt.legend)[i])
			i += 1

	## Choose a filename and save the Matlab script
	def saveAs(self):
		fileName = QtGui.QFileDialog.getSaveFileName(self)
		self.lineEdit_4.clear()
		self.lineEdit_4.insert(fileName)
		self.save()

	## Save the Matlab script
	def save(self):
		if self.tt==0:
			return
		fileName = self.lineEdit_4.text()
		self.tt.Write_Matlab_Script(fileName)
		print "Matlab file succesfully created."

	## Exit the program
	def exit(self):
		configs = open ('config.txt', 'r')
		f = configs.readlines()
		f[1] = str(self.lineEdit_db.text()) + "\n"
		f[3] = str(self.lineEdit_6.text()) + "\n"
		f[5] = self.lineEdit_4.text() + "\n"
		f[7] = str(self.comboBox_4.currentIndex()) + '\n'
		f[9] = str(self.comboBox_5.currentIndex()) + "\n"
		self.writeConfigFile(f)
		if not (self.tt==0):
			self.tt.Close()
			del self.tt
		sys.exit()

	## Write GUI config file (config.txt)
	def writeConfigFile(self, f):
		configs = open ('config.txt', 'w')
		y = len(f)-1
		i = 0
		while i <= y:
			configs.write(f[i])
			i += 1
		configs.close()
		print "Configuration saved succesful."

	## Ask if the GUI shall be closed
	def closeEvent(self, event):
		quit_msg = "Are you sure you want to exit the GUI?"
		reply = QtGui.QMessageBox.question(self, 'Message', 
		quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
		if reply == QtGui.QMessageBox.Yes:
			self.exit()
			event.accept()
		else:
			event.ignore()

## Application
app = QtGui.QApplication(sys.argv) 
## Layout
dialog = Gui() 
dialog.show() 
sys.exit(app.exec_())


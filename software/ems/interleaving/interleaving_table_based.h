//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * interleaving_table_based.h
 *
 *  Created on: 15.02.2011
 *      Author: gimmler
 */

#ifndef INTERLEAVING_TABLE_BASED_H_
#define INTERLEAVING_TABLE_BASED_H_

#include <fstream>
#include "../modules/interleaving/interleaving_iface.h"
#include "interleaving_table_based_param.h"
#include "interleaving_table_based.h"

using namespace cse_lib;

/**
 * \brief simple interleaver which reads index from a file
 *
 * This is a very basic interleaver class. It needs a given interleaver index in a file.
 * First, it will try to open a stored interleaver file given by the param interleaver_table.
 * If NO_TABLE is chosen, it will try to open the file given by the param file_name.
 * If that file does not exist it generates the index "0, 1, 2, 3, ..." etc. for the necessary length.
 * Text files containing interleaver tables start with the length of the table followed by all entries.
 */
template<class T>
class Interleaving_Table_Based: public Interleaving_Table_Based_Parameter, public Interleaving_Interface<T> {
public:
	Interleaving_Table_Based();
	virtual ~Interleaving_Table_Based();

	///process one block of input data
	int Run();

protected:
	///read interleaver table and resize memories
	void Init();
	///resize all memories
	void Alloc_Mem();
	///open the chosen interleaver table text file
	void open_interleaver_table_file(ifstream *filepointer);

private:
	///stored interleaver table
	Buffer<int> interleaver_index;
	///length of interleaver table
	int interleaver_length;

};

#endif /* INTERLEAVING_TABLE_BASED_H_ */

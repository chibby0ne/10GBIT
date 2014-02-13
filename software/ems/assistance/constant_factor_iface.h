//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * constant_factor_iface.h
 *
 *  Created on: 09.01.2012
 *      Author: gimmler
 */

#ifndef CONSTANT_FACTOR_IFACE_H_
#define CONSTANT_FACTOR_IFACE_H_

#include "../modules/base/base_iface.h"

using namespace cse_lib;

/// Interface class for puncturing units
/**
 * \ingroup interfaces
 */
template<class T> class Constant_Factor_Interface: public Base_Interface {
public:
	Constant_Factor_Interface() {
		input.Register("input", input_data_list_);
		output.Register("output", output_data_list_);
	}

	virtual ~Constant_Factor_Interface() {
	}

	///input buffer
	Data_In<T> input;

	///output buffer = factor * input buffer
	Data_Out<T> output;

};

// Explicit instantiation
template class Constant_Factor_Interface<int> ;
template class Constant_Factor_Interface<float> ;
template class Constant_Factor_Interface<double> ;
template class Constant_Factor_Interface<complex<int> >;
template class Constant_Factor_Interface<complex<float> >;
template class Constant_Factor_Interface<complex<double> >;

#endif /* CONSTANT_FACTOR_IFACE_H_ */

//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * mimo_channel_rayleigh_iface.h
 *
 *  Created on: 03.02.2011
 *      Author: gimmler
 */

#ifndef MIMO_CHANNEL_RAYLEIGH_IFACE_H_
#define MIMO_CHANNEL_RAYLEIGH_IFACE_H_

#include "../modules/base/base_iface.h"

using namespace cse_lib;

///interface class for Rayleigh fading mimo channel
template<class T>
class Mimo_Channel_Rayleigh_Interface: public Base_Interface {
public:
	Mimo_Channel_Rayleigh_Interface() {
		input_symb.Register("input_symbols", input_data_list_);
		output_symb.Register("output_symbols", output_data_list_);
		output_channel_matrices.Register("output_channel_matrices", output_data_list_);
	}
	virtual ~Mimo_Channel_Rayleigh_Interface() {
	}

	///input modulation symbols
	Data_In<complex<T> > input_symb;
	///output received symbols
	Data_Out<complex<T> > output_symb;
	///output all generated channel matrices of one block
	Data_Out<complex<T> , 2> output_channel_matrices;

};

#endif /* MIMO_CHANNEL_RAYLEIGH_IFACE_H_ */

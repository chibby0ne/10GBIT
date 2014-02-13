//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * dec_conv_iface.h
 *
 *  Created on: 16.02.2011
 *      Author: gimmler
 */

#ifndef DEC_CONV_IFACE_H_
#define DEC_CONV_IFACE_H_

#include "../modules/base/base_iface.h"

///interface for NSC convolutional decoder
template<class T>
class Decoder_Convolutional_Interface: public Base_Interface {

public:

	Decoder_Convolutional_Interface() {
		input_llrs_apriori.Register("input_llrs_apriori", input_data_list_);
		output_bits.Register("output_bits", output_data_list_);
		output_llrs_extrinsic.Register("output_llrs_extrinsic", output_data_list_);
		output_llrs_app.Register("output_llr_app", output_data_list_);
		valid_codeword_differences.Register("valid_codeword_differences", status_out_list_);
	}

	virtual ~Decoder_Convolutional_Interface() {
	}

	///input llrs apriori information of the whole codeword
	Data_In<T> input_llrs_apriori;
	///output llrs extrinsic information of the whole codeword
	Data_Out<T> output_llrs_extrinsic;
	///output llrs a-posteriori information of the whole codeword
	Data_Out<T> output_llrs_app;
	///output bits of the systematic bits
	Data_Out<unsigned int, 2> output_bits;
	///number of differences for re-encoding, if 0 a valid codeword has been found
	Status_Out<unsigned int> valid_codeword_differences;
};

#endif /* DEC_CONV_IFACE_H_ */

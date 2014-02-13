//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * enc_conv.h
 *
 *  Created on: 15.02.2011
 *      Author: gimmler
 */

#ifndef ENC_CONV_H_
#define ENC_CONV_H_

#include "enc_conv_param.h"
#include "conv_share.h"
#include "../modules/encoder/enc_iface.h"

///convolutional NSC encoder with termination in zero-state
class Encoder_Convolutional: public Encoder_Interface, public Encoder_Convolutional_Parameter, public Convolutional_Share {
public:
	Encoder_Convolutional() {
	}
	virtual ~Encoder_Convolutional() {
	}
	///process one block of input data
	int Run();

protected:
	///compute derived parameters and resize memories
	void Init();
	///resize all memories
	void Alloc_Mem();
	///returns the length of one code word for the current polynomials and the infoword length
	unsigned int get_coded_bit_length(unsigned int uncoded_bit_length);

private:
	unsigned int uncoded_bit_length_;
};

#endif /* ENC_CONV_H_ */

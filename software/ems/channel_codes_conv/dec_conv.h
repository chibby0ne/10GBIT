//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * dec_conv.h
 *
 *  Created on: 16.02.2011
 *      Author: gimmler
 */

#ifndef DEC_CONV_H_
#define DEC_CONV_H_

#include "dec_conv_param.h"
#include "conv_share.h"
#include "dec_conv_iface.h"

///convolutional NSC decoder with termination in zero-state
/**
 * BCJR / MAP decoding of non-systematic, non-recursive convolutional codes.
 */
template<class T>
class Decoder_Convolutional: public Decoder_Convolutional_Interface<T> , public Decoder_Convolutional_Parameter, public Convolutional_Share {

public:
	Decoder_Convolutional();
	virtual ~Decoder_Convolutional();
	///process one block of input data
	int Run(int Iteration);

	///process one block of input data
	int Run();

protected:
	///compute derived parameters and resize memories
	void Init();
	///resize all memories
	void Alloc_Mem();
	///returns the length of the trellis for the current polynomials and the codeword length
	unsigned int get_trellis_length(unsigned int coded_bit_length);

	///determines Gammas for a set of parity information
	void calc_Gammas(T par1, T par2, Mat<T> & Gammas);
	///bcjr algorithm
	void bcjr_algorithm(int Iteration = 0);
	///internal function which clips all extrinsic output LLR values
	void ClipLlrs();
	///check if output of decode is a valid codeword, returns 0 if codeword is valid, otherwise number of differences
	int check_codeword();
	///max and max_star operation
	inline T max_star(T a, T b);

private:
	///length of one codeword
	unsigned int coded_bit_length_;
	///length of the trellis = number of sytematic bits + termination bits
	int trellis_length_;
};

#endif /* DEC_CONV_H_ */

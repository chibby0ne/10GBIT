//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the demapper class.
/// \author Timo Lehnigk-Emden
/// \date   2010/12/10
//
#ifndef DEMAPPER_H_
#define DEMAPPER_H_

#include <complex>
#include <iostream>
#include "demapper_iface.h"
#include "demapper_param.h"
#include "mapper_share.h"


namespace cse_lib {

/// SISO Signal Demapper class
/**
 * This module demaps the received symbols from the channel and extracts the
 * hard bits, the LLR values, and the bit probabilities.
 *
 * \ingroup modules
 */
class Demapper: public Demapper_Interface,
                public Demapper_Parameter,
                public Mapper_Share
{

public:

	Demapper();
	virtual ~Demapper();

	/// Demaps one block of data.
	int Run();

private:

	void Alloc_Mem();

	void Init();

	void Demapping_BPSK(Buffer<complex<float> >& input_symb, Buffer<float>& output_bits_llr);
	void Demapping_QPSK(Buffer<complex<float> >& input_symb, Buffer<float>& output_bits_llr);

	void Sub_Optimal_Demapping_16QAM(Buffer<complex <float> >& input_symb,
	                                 Buffer<float>& output_bits_llr);

	void Sub_Optimal_Demapping_64QAM(Buffer<complex <float> >& input_symb,
	                                 Buffer<float>& output_bits_llr);

	void Optimal_Demapping(Buffer<complex<float> >& input_symb,
	                       Buffer<float, 2>& output_bits_prob,
	                       Buffer<float>& output_bits_llr);

	void Demapping_8PSK(Buffer<complex<float> >& input_symb,
	                    Buffer<float>& output_bits_llr);

	void Init_Demapping();


	/// Number of symbols per constellation point.
	unsigned int num_of_symbols_;

	unsigned int bits_per_symbol_;

	Buffer<int,2> S0_;
	Buffer<int,2> S1_;

	void Dump_Wrong_Constellation_Message();

	void Calc_Hard_Bits();
	void Scale_CRF(Buffer<float> &output_bits_llr);
};
}
#endif


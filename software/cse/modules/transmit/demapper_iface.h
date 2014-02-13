//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the data interface class of the SISO demapper.
/// \author Timo Lehnigk-Emden
/// \date   2010/12/10
//
#ifndef DEMAPPER_IFACE_H_
#define DEMAPPER_IFACE_H_

#include <complex>
#include "../base/base_iface.h"


namespace cse_lib {

/// SISO demapper interface class
/**
 * This class defines the interface for SISO demappers. The input port expects
 * the symbols from a channel. The results are the demapped bits (hard, LLR,
 * probabiliety values).
 *
 *  \ingroup interface
 */
class Demapper_Interface: public Base_Interface
{

public:

	Demapper_Interface()
	{
		input_symb.Register("input_symb", input_data_list_);
		output_bits.Register("output_bits", output_data_list_);
		output_bits_llr.Register("output_bits_llr", output_data_list_);
		output_bits_prob.Register("output_bits_prob", output_data_list_);

		histogram.Register("histogram", status_out_list_);
	}

	virtual ~Demapper_Interface(){}

	/// Input symbol port, for the received symbols
	Data_In<complex<float> > input_symb;

	/// Output hard bits port, demapped bits
	Data_Out<unsigned int> output_bits;

	/// Output LLR port, LLR values of the demapped bits
	Data_Out<float> output_bits_llr;

	/// Output probability port, probabilities of the demapped bits
	Data_Out<float, 2> output_bits_prob;

	/// Histrogram showing the saturation degree
	Status_Out<float, 0, Status_Out_Plugin_Histogram> histogram;
};
}
#endif

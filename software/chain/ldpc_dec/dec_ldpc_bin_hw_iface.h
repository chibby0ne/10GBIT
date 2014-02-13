//
//  Copyright (C) 2010 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Channel decoder interface
/// \author Matthias Alles
/// \date   2010/12/03
//

#ifndef DEC_LDPC_BIN_HW_IFACE_H_
#define DEC_LDPC_BIN_HW_IFACE_H_

#include "../../cse/modules/base/base_iface.h"

namespace cse_lib {

/// Binary LDPC Decoder interface class
/**
 * \ingroup interface
 */
template<typename T>
class Decoder_LDPC_Binary_HW_Interface: public Base_Interface
{

public:

	Decoder_LDPC_Binary_HW_Interface()
	{
		input_bits_llr.Register("input_bits_llr", input_data_list_);
		output_bits.Register("output_bits", output_data_list_);
		output_bits_llr_app.Register("output_bits_llr_app", output_data_list_);

		mean_iterations.Register("mean_iterations", status_out_list(), true);
		iterations_performed.Register("iterations_performed", status_out_list(), false);
		flipped_bits.Register("flipped_bits", status_out_list_, true);
		decoding_successful.Register("decoding_successful", status_out_list(), false);
		num_modified_systematic_bits.Register("num_modified_systematic_bits", status_out_list(), false);
	}

	virtual ~Decoder_LDPC_Binary_HW_Interface() {}

	/// Quantized channel values.
	Data_In<T>               input_bits_llr;

	/// Hard bits after decoding, without shortened and parity bits. Not truncated, information bits have to be extraced by hand.
	Data_Out<unsigned int, 2> output_bits;

	/// LLR values of decoded bits.
	Data_Out<T, 2>            output_bits_llr_app;

	/// Amount of iterations required for decoding.
	Status_Out<unsigned int>  iterations_performed;

	/// Number of flipped bits per iterations (for each single iteration number).
    Status_Out<unsigned int, 1> flipped_bits;

	/// Mean numer of iterations required for decoding (for each single iteration number).
	Status_Out<unsigned int, 1, Status_Out_Plugin_Mean> mean_iterations;

	/// All parity checks satisfied?
	Status_Out<bool>          decoding_successful;

	/// Number of modified systematic bits while decoding.
	Status_Out<unsigned int>  num_modified_systematic_bits;

	/// Number of parity checks that are not satisfied after decoding.
	Status_Out<unsigned int>  num_unsatisfied_parity_checks;
};
}
#endif // DEC_LDPC_BIN_HW_IFACE_H_

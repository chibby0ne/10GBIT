//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  This is the interface definition for the error rates statistics class
/// \author Timo Lehnigk-Emden
/// \date   2010/12/22
//

#ifndef STAT_ERROR_RATES_IFACE_H_
#define STAT_ERROR_RATES_IFACE_H_

#include "../base/base_iface.h"

namespace cse_lib {

/// Data interface definition for the error rates statistics class
/**
 * \ingroup interface
 */
template <unsigned int DIM> class Statistics_Error_Rates_Interface : public Base_Interface
{

public:

	Statistics_Error_Rates_Interface()
	{
		input_bits.Register("input_bits", input_data_list_);
		input_bits_ref.Register("input_bits_ref", input_data_list_);

		num_total_bits.Register("num_total_bits", status_out_list_);
		num_diff_bits.Register("num_diff_bits", status_out_list_);
		error_rate_bits.Register("error_rate_bits", status_out_list_);
		num_total_blocks.Register("num_total_blocks", status_out_list_);
		num_diff_blocks.Register("num_diff_blocks", status_out_list_);
		error_rate_blocks.Register("error_rate_blocks", status_out_list_);
	}

	virtual ~Statistics_Error_Rates_Interface() {}

	/// Input bits or symbols
	Data_In<u_int, DIM> input_bits;

	/// Reference bits or symbols
	Data_In<u_int> input_bits_ref;

	/// Total number of bits read
	Status_Out<ull_int> num_total_bits;

	/// Total number of different bits
	Status_Out<ull_int, DIM-1, Status_Out_Plugin_Default> num_diff_bits;

	/// Error rate of the bits
	Status_Out<float, DIM-1, Status_Out_Plugin_Default> error_rate_bits;

	/// Total number of blocks or frames read
	Status_Out<ull_int> num_total_blocks;

	/// Total number of different blocks or frames
	Status_Out<ull_int, DIM-1, Status_Out_Plugin_Default> num_diff_blocks;

	/// Error rate of the blocks or frames
	Status_Out<float, DIM-1, Status_Out_Plugin_Default> error_rate_blocks;
};
}
#endif

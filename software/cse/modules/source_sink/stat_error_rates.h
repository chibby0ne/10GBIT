//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of error rates statistics class
/// \author Timo Lehnigk-Emden
/// \date   2010/12/21
//

#ifndef STAT_ERROR_RATES_H_
#define STAT_ERROR_RATES_H_

#include "stat_error_rates_iface.h"
#include "stat_error_rates_param.h"
#include "stat_error_rates_share.h"

namespace cse_lib
{

/// Error rates statistics class and simulation termination (sink)
/**
 * \ingroup modules
 */
template <unsigned int DIM>
class Statistics_Error_Rates : public Statistics_Error_Rates_Interface<DIM>,
                               public Statistics_Error_Rates_Parameter,
                               public Statistics_Error_Rates_Share
{

public:

	Statistics_Error_Rates() {};
	virtual ~Statistics_Error_Rates() {};

	/// Process one block
	int Run();

	/// Print the current results of the error rate calculation
	//void Dump_Intermediate_Result();

	// Make inherited parameter from template class visible
	using Statistics_Error_Rates_Interface<DIM>::input_bits;
	using Statistics_Error_Rates_Interface<DIM>::input_bits_ref;
	using Statistics_Error_Rates_Interface<DIM>::num_diff_bits;
	using Statistics_Error_Rates_Interface<DIM>::error_rate_bits;
	using Statistics_Error_Rates_Interface<DIM>::num_total_blocks;
	using Statistics_Error_Rates_Interface<DIM>::num_total_bits;
	using Statistics_Error_Rates_Interface<DIM>::num_diff_blocks;
	using Statistics_Error_Rates_Interface<DIM>::error_rate_blocks;
	using Statistics_Error_Rates_Interface<DIM>::input_data_list_;

private:

	void Init();

	Buffer<ull_int> num_total_blocks_;
	Buffer<ull_int> num_total_bits_;
	Buffer<ull_int> num_diff_blocks_;
	Buffer<ull_int> num_diff_bits_;

	unsigned int num_iter_;
};


/// Specialization for one dimensional input
template <>
class Statistics_Error_Rates<1> : public Statistics_Error_Rates_Parameter,
                                  public Statistics_Error_Rates_Interface<1>,
                                  public Statistics_Error_Rates_Share
{

public:

	Statistics_Error_Rates() {};
	virtual ~Statistics_Error_Rates() {};

	/// Process one block
	int Run();

	// Make inherited parameter from template class visible
	using Statistics_Error_Rates_Interface<1>::input_bits;
	using Statistics_Error_Rates_Interface<1>::input_bits_ref;
	using Statistics_Error_Rates_Interface<1>::num_diff_bits;
	using Statistics_Error_Rates_Interface<1>::error_rate_bits;
	using Statistics_Error_Rates_Interface<1>::num_total_blocks;
	using Statistics_Error_Rates_Interface<1>::num_total_bits;
	using Statistics_Error_Rates_Interface<1>::num_diff_blocks;
	using Statistics_Error_Rates_Interface<1>::error_rate_blocks;
	using Statistics_Error_Rates_Interface<1>::input_data_list_;

private:

	void Init();

	ull_int num_total_blocks_;
	ull_int num_total_bits_;
	ull_int num_diff_blocks_;
	ull_int num_diff_bits_;
};

}
#endif

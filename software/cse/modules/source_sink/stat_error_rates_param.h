//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  This is the parameter definition for the error rates statistics class
/// \author Timo Lehnigk-Emden
/// \date   2010/12/21
//

#ifndef STAT_ERROR_RATES_PARAM_H_
#define STAT_ERROR_RATES_PARAM_H_

#include "../base/base_param.h"

namespace cse_lib {

/// Parameter definition for the error rates statistics class
/**
 * \ingroup parameter
 */
class Statistics_Error_Rates_Parameter : public Base_Parameter
{

public:

	Statistics_Error_Rates_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Statistics_Error_Rates_Parameter() {}

	/// Unique identification string of the module
	static string Unique_ID(){ return "Statistics_Error_Rates";}


	/**************
	 * Parameters *
	 **************/

	/// Number of total blocks to process
	/**
	 * If this number is reached the module's Run function returns a value to
	 * indicate that the maximum number is reached.
	 */
	Param<ull_int> max_num_total_blocks;

	/// Number of total blocks which are different
	/**
	 * If this number is reached the module's run function returns a value to
	 * indicate that the maximum number is reached.
	 */
	Param<ull_int> max_num_diff_blocks;

	/// Controls the console output behavior
	/**
	 * Print the current error rates after each block on the console if true (=default).
	 */
	Param<bool> print_status_permanent;

	/// Name of inner port dimensions
	/**
	 * Set the name of inner dimension of the status out ports.
	 * This name is important for the results file.
	 * If you want to track the error rates over different iterations for example.
	 * The name should be "iteration" that is the default value.
	 */
	Param<raw_string> out_port_inner_dim_name;

	/// Output address offset
	/**
	 * This value is the offset between the input buffer dimension and the status out port results
	 * dimension. If for example the dimension "0" of the input buffer contains the data of the first
	 * iterations. Then the value should be "1" which is the default value, to ensure that in
	 * result file the result data have the dimension number "1" and not "0".
	 */
	Param<u_int> out_port_inner_dim_addr_offset;

protected:


	/******************
	 * Default values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - max_num_diff_blocks  : 100
	 *  - max_num_total_blocks : 1000
	 */
	void Set_Default_Values()
	{
		max_num_diff_blocks.Init(100, "max_num_diff_blocks", param_list_);
		max_num_total_blocks.Init(1000, "max_num_total_blocks", param_list_);
		print_status_permanent.Init(true, "print_status_permanent", param_list_);
		out_port_inner_dim_name.Init("iteration", "out_port_inner_dim_name", param_list_);
		out_port_inner_dim_addr_offset.Init(1, "out_port_inner_dim_addr_offset", param_list_);
	}
};
}
#endif

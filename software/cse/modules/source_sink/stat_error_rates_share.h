//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Shared functionality for the error rates counter class.
/// \author Timo Lehnigk-Emden
/// \date   2010/12/22
//

#ifndef STAT_ERROR_RATES_SHARE_H_
#define STAT_ERROR_RATES_SHARE_H_

#include "../base/cse_namespace.h"
#include <iomanip>
#include <ctime>

namespace cse_lib {


/// Shared functionalities for the statistics error rates class.
/**
 * \ingroup share
 */
class Statistics_Error_Rates_Share
{
public:

	Statistics_Error_Rates_Share() {}
	virtual ~Statistics_Error_Rates_Share() {}

protected:

	/// Count the number of differences (errors).
	/**
	 * Compares the two buffers and updates the given parameters. In case the
	 * two buffers differ in length, the shortest length is taken as reference.
	 *
	 * \param data Input data.
	 * \param ref  Reference data.
	 * \param num_total_blocks variable to store the number of total blocks processed
	 * \param num_total_bits variable to store the number of total bits processed
	 * \param num_diff_blocks variable to store the number of blocks differs
	 * \param num_diff_bits variable to store the number of total bits differs
	 */
	void Count_Errors(Buffer<u_int>& data, Buffer<u_int>& ref,
	                  ull_int& num_total_blocks, ull_int& num_total_bits,
	                  ull_int& num_diff_blocks, ull_int& num_diff_bits);


	/// Print the current status of the error counter unit
	void Dump_Intermediate_Result(string instance_name,
	                              ull_int max_num_diff_blocks, ull_int max_num_total_blocks,
	                              ull_int num_diff_blocks, ull_int num_total_blocks,
	                              float error_rate_blocks, float error_rate_bits, float thrpt = -1.0)
	{
		cout << "\r" << instance_name << ": " << setw(10) << right << num_total_blocks << " frames read - ("
		     << static_cast<int> (100 * (num_total_blocks /
		        static_cast<float> (max_num_total_blocks)))
		     << "%)" << "("
		     << static_cast<int> (100 * (num_diff_blocks /
		        static_cast<float> (max_num_diff_blocks)))
		     << "%)"
		     << " - FER: " << setprecision(2) << scientific << error_rate_blocks
		     << " - BER: " << setprecision(2) << scientific << error_rate_bits;

		cout.unsetf(cout.scientific);

		if(thrpt > 0)
			cout << " - Speed: " << fixed <<setprecision(2) << thrpt << " kbit/s";

		cout << flush;
	}

	/// Calculates the throughput between two consecutive calls of the function.
	/**
	 * Before utilization, set init_time_ to clock(), since it initiates the
	 * start point of measurement.
	 */
	float Calc_Throughput(ull_int amount_of_data)
	{
		float ret;
		clock_t curr_time = clock();
		static float tmp_val = CLOCKS_PER_SEC / 1000;

		clock_t diff = curr_time - init_time_;

		if(diff > 0)
			ret =  static_cast<float> (amount_of_data) * tmp_val / diff;
		else
			ret =  0.0;

		return ret;
	}

	/// Set this time to clock() when configuring in order to use Calc_Throughput() function.
	clock_t init_time_;

};

}

#endif

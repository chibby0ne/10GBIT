//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Mean plug-in for the status out port
/// \author Timo Lehnigk-Emden
/// \date   2011/02/04
//

#ifndef STATUS_OUT_PLUGIN_MEAN_H_
#define STATUS_OUT_PLUGIN_MEAN_H_

#include "status_out_plugin_abstract.h"

namespace cse_lib{
/// Mean value plug-in for the Status_Out port.
/**
 * This class calculates the mean value of all values written.
 * For this an internal sum is built of all values written to one location.
 * Furthermore, the number of write operations is stored.
 *
 * In case of a Read(), the quotient, i.e., the mean value is calculated.
 *
 * \ingroup status_out_plugin
 */
template<class T>
class Status_Out_Plugin_Mean : public Status_Out_Plugin_Abstract
{
public:

	Status_Out_Plugin_Mean()
	{
		Reset();
		valid_flag_ = true;
	}

	/// Output data type of the plugin.
	typedef float OUTPUT_TYPE;

	/// Reset the plug-in
	void Reset()
	{
		sum_ = 0;
		number_ = 0;
	}

	/// Return the mean value
	OUTPUT_TYPE Read()
	{
		return static_cast<OUTPUT_TYPE> (sum_) / number_;
	}

	/// Write the value
	/**
	 *  \param value to add the the mean value
	 */
	void Write(T value)
	{
		sum_ += value;
		number_++;
	}

	/// Gives one datum of the internal Buffer and returns true as long as data is available.
	bool Get_Value_Str(Status_Out_Data_Struct& data)
	{
		if(number_ != 0 && valid_flag_ == true)
		{
			data.value = Conv_Num_To_Str(Read());
			valid_flag_ = false;
			return true;
		}
		else
		{
			valid_flag_ = true;
			return false;
		}
	}

private:
	// Sum memory
	T sum_;
	// Number of input values
	ull_int number_;

	bool valid_flag_;
};
}

#endif

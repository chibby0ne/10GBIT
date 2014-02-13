//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Histogramm plug-in for the statuts out port
/// \author Timo Lehnigk-Emden
/// \date   2011/02/04
//

#ifndef STATUS_OUT_PLUGIN_HIST_H_
#define STATUS_OUT_PLUGIN_HIST_H_

#include "status_out_plugin_abstract.h"

namespace cse_lib {
/// Histogram plug-in for the Status_Out port.
/**
 * This plug-in counts, how many times a certain value has occurred.
 * The generated histogram is normalized to the overall number of values, such
 * that the output values for every single bin are always between 0.0 and 1.0.
 * The address for the value is calculated by value / bin_size. By default,
 * the bin size is 1.
 *
 * \ingroup status_out_plugin
 */
template<class T>
class Status_Out_Plugin_Histogram : public Status_Out_Plugin_Abstract
{
public:
	Status_Out_Plugin_Histogram()
	{
		Reset();
	}

	virtual ~Status_Out_Plugin_Histogram() {}

	/// Data type of the output values
	typedef float OUTPUT_TYPE;

	/// Reset the memory.
	void Reset()
	{
		curr_pos_ = 0;
		dim_name_ = "bin";
		counter_ = 0;
		bin_size_ = 1;
		hist_positive_.Resize(1);
		hist_negative_.Resize(1);
		hist_positive_.Clear();
		hist_negative_.Clear();
	}

	/// Return the histogram value at the position index.
	OUTPUT_TYPE Read(int index)
	{
		if (index < 0)
		{
			unsigned int buf_addr = static_cast<unsigned int> (index * -1);
			if (hist_negative_.length() > 0)
				return static_cast<OUTPUT_TYPE> (hist_negative_[buf_addr]) / counter_;
			else
				return 0.0;
		}
		else
		{
			unsigned int buf_addr = static_cast<unsigned int> (index);
			if (hist_positive_.length() > 0)
				return static_cast<OUTPUT_TYPE> (hist_positive_[buf_addr]) / counter_;
			else
				return 0.0;
		}
	}

	/// Increment the histogram at position value.
	/**
	 * If the position is out of range the memory will be automatically
	 * extended and the previous values will be kept.
	 */
	void Write(T value)
	{
		int index = static_cast<int> (value / bin_size_);
		if (index >= 0)
		{
			unsigned int buf_addr = static_cast<unsigned int> (index);

			if (hist_positive_.length() <= buf_addr)
				hist_positive_.Resize(buf_addr + 1, true);
			hist_positive_[buf_addr]++;
		}
		else
		{
			unsigned int buf_addr = static_cast<unsigned int> (index * -1);

			if (hist_negative_.length() <= buf_addr)
				hist_negative_.Resize(buf_addr + 1, true);
			hist_negative_[buf_addr]++;
		}

		counter_++;
	}

	/// Return the number of histogram points which are currently stored.
	u_int length()
	{
		// The negative histogram holds one value for the zero as well, which is not used.
		return hist_positive_.length() + hist_negative_.length() - 1;
	}

	/// Gives one datum of the internal Buffer and returns true as long as data is available.
	bool Get_Value_Str(Status_Out_Data_Struct& data)
	{
		if(curr_pos_ < length())
		{
			// Get current position as negative address.
			int curr_pos_int = static_cast<int>(curr_pos_ - hist_negative_.length() + 1);

			// Empty data structure before writing.
			data.Clear();

			data.value = Conv_Num_To_Str(Read(curr_pos_int));
			string str_tmp = Conv_Num_To_Str(curr_pos_int);
			data.address.Append(str_tmp);
			data.dim_name.Append(dim_name_);

			curr_pos_++;
			return true;
		}
		else
		{
			curr_pos_ = 0;
			return false;
		}
	}

	/// Set the distance between two bins of the histogram, default is 1.
	void bin_size(T bin_size)
	{
		bin_size_ = bin_size;
	}

private:

	unsigned int curr_pos_;
	string dim_name_;

	ull_int counter_;
	T bin_size_;

	// Negative histogram, has also place for a zero, which is not used.
	Buffer<ull_int> hist_negative_;
	Buffer<ull_int> hist_positive_;
};
}
#endif

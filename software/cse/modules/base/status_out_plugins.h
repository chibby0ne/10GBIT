//
//  Copyright © 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems
//
/// \file
/// \brief  Plug-in interface and default plug-in for the status out port
/// \author Timo Lehnigk-Emden
/// \date   2011/01/11
//

#ifndef STATUS_OUT_PLUGINS_H_
#define STATUS_OUT_PLUGINS_H_

#include "status_out_plugin_abstract.h"

namespace cse_lib {

/// Default plug-in for the Status_Out port (memory behaviour).
/**
 * This plug-in is used by Status_Out if no other plug-in is selected.
 * The default plug-in works as a simple memory for one value.
 *
 * Each Write() call overwrites the current value.
 * The Read() function returns the currently stored value.
 *
 * \ingroup status_out_plugin
 */
template<class T>
class Status_Out_Plugin_Default : public Status_Out_Plugin_Abstract
{
public:

	Status_Out_Plugin_Default()
	{
		Reset();
	}

	/// Reset the memory
	void Reset()
	{
		value_ = 0;
		valid_flag_ = true;
	}

	/// Read the stored value
	T Read()
	{
		return value_;
	}

	/// Store one value in the plug-in
	void Write(T value)
	{
		value_ = value;
	}

	/// Write the value in the given data structure
	/**
	 * \return true: data is valid, false: data is invalid
	 */
	bool Get_Value_Str(Status_Out_Data_Struct& data)
	{
		if(valid_flag_ == true)
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

	T value_;
	bool valid_flag_;
};
}
#endif

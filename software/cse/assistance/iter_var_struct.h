//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Date structure for iteration variables in the configuration
/// \author Timo Lehnigk-Emden
/// \date   2010/12/20
//

#ifndef ITER_VAR_STRUCT_H_
#define ITER_VAR_STRUCT_H_

#include <string>
#include <sstream>
#include "generic_list.h"

/// Data structure to hold the values of an iteration variable of the configuration.
/**
 * This class is a data structure to hold iterative values given by a string.
 * All values are stored in a list of strings.
 * The class provides a convenient way to set the iteration space with the function Set_From_String()
 * \ingroup base_class
 */
class Iteration_Variable_Structure
{

public:

	Iteration_Variable_Structure()
	{
		curr_val_.clear();
		curr_index_ = 0;
	}

	/// Create a new instance and initialize the values with the given string and the name
	/**
	 * \param iter_str string with the iterative value (details see Set_From_String())
	 * \param var_name name of the instance
	 */
	Iteration_Variable_Structure(std::string iter_str, std::string var_name)
	{
		Set_From_String(iter_str, var_name);
	}

	virtual ~Iteration_Variable_Structure() { }


	/// Set the value from a given string and set the instance name
	/**
	 * \param iter_str string which contains the iteration string
	 * The values can be semicolon separated, a range can be given as "start:incr:end" or "start:end".
	 * For the last case an increment of one is assumed.
	 * \param var_name the name to set of the instance (name of the iteration variable)
	 */
	void Set_From_String(std::string iter_str, std::string var_name)
	{
		iter_var_list_.Clear();
		name_ = var_name;
		//
		// String processing
		//
		std::string result_temp;
		size_t pos, pos3;

		do
		{
			// Find separator for the enumeration of different values
			pos = iter_str.find(';');

			// Get String part before the comma without the comma
			result_temp = iter_str.substr(0, pos);

			// Get String part after the comma without the comma
			iter_str = iter_str.substr(pos + 1);

			// Detect : for range description start:incr:end or start:end with incr = 1
			pos3 = result_temp.find(':');

			if (pos3 != std::string::npos)
			{
				double start, end, incr = 1.0;

				// Get start value
				hlp_fct::str::Convert_Strict(result_temp.substr(0, pos3), start);

				size_t pos_tmp;

				result_temp = result_temp.substr(pos3 + 1);

				// Search for the second :
				pos_tmp = result_temp.find(':');

				// Found it
				if (pos_tmp != std::string::npos)
				{
					// Extract increment value
					hlp_fct::str::Convert_Strict(result_temp.substr(0, pos_tmp), incr);

					// Update string
					result_temp = result_temp.substr(pos_tmp + 1);
				}

				// Read end value
				hlp_fct::str::Convert_Strict(result_temp, end);

				// Check for values which results in infinity number of points
				Check_Consistency(var_name, start, incr, end);

				// Iterate over the range with the increment
				if(incr != 0.0)
				{
					double v=start;
					std::string tmp_str;

					// Handle the value precision problem with an incremental value.
					for(unsigned int i = 0; v <= end; i++, v = start + i*incr)
					{
						tmp_str = (hlp_fct::str::Conv_Num_To_Str(v));
						iter_var_list_.Append(tmp_str);
					}
				}
				else
				{
					std::string tmp_str(hlp_fct::str::Conv_Num_To_Str(start));
					iter_var_list_.Append(tmp_str);
				}
			}

			// Read single value if no : character is present in the string
			if((pos3 == std::string::npos))
			{
				//TODO: Add exception in case of empty string?
				iter_var_list_.Append(result_temp);
			}
		}
		while (pos != std::string::npos);

		curr_index_ = 0;
		if(iter_var_list_.length() > 0)
			curr_val_ = iter_var_list_[0];
		else
			curr_val_.clear();
	}

	/// Update the values to the next iteration
	int Update_To_Next_Value()
	{

		// Last value reached, range check
		if(curr_index_ == iter_var_list_.length() - 1)
		{
			curr_index_ = 0;
			curr_val_ = iter_var_list_[curr_index_];
			return -1;
		}

		curr_index_++;
		curr_val_ = iter_var_list_[curr_index_];
		return 0;
	}

	/// Check value consistency
	void Check_Consistency(std::string var_name, double& start_val, double& incr_val, double& end_val)
	{
		// Check for zero incremental
		if( incr_val == 0.0 && end_val != start_val)
		{
			hlp_fct::logging::Msg(hlp_fct::logging::WARNING, "Iteration_Variable_Structure",
			                    "Increment value of variable with name \"" + var_name + "\" "
			                    "is zero and start value, end value are different. "
			                    "Setting end value to start value to prevent infinity runtime!");
			end_val = start_val;
		}

		// Check for unreachable end point (diverged parameter)
		if(fabs(start_val + incr_val - end_val) > fabs(start_val - end_val))
		{
			hlp_fct::logging::Msg(hlp_fct::logging::WARNING, "Iteration_Variable_Structure",
			                    "Parameter ensemble of variable with name \"" + var_name + "\" "
			                    "leads to an infinity runtime!. Only the start value is taken "
			                    "and the increment value is set to zero!");
			end_val = start_val;
			incr_val = 0.0;
		}
	}

	/// Return the current value of the list
	std::string curr_val() { return curr_val_;}

	/// Return the name of the instance
	std::string& name() { return name_;}

	/// Set the name of the instance
	void name(std::string name) {name_ = name;}

	/// Return the number of iteration points
	unsigned int length() {return iter_var_list_.length();}

	/// Streaming operator to print out the values
	/**
	 * Elements are separated by semicolon. That format corresponds to input format of the
	 * Set_From_String data input function.
	 */
	inline friend std::ostream & operator<<(std::ostream& os, Iteration_Variable_Structure& data)
	{
		//os << data.iter_var_list_;
		unsigned int i;
		for(i = 0; i < data.iter_var_list_.length() - 1; i++)
			os << data.iter_var_list_.Get_Ref(i) << "; ";

		os << data.iter_var_list_.Get_Ref(i);

		return os;
	}

private:

	// List to store all data points
	Generic_List<std::string, std::vector> iter_var_list_;

	// current value
	std::string curr_val_;

	// current index
	unsigned int curr_index_;

	// name of the variable
	std::string name_;
};
#endif

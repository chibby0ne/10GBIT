//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Parameter template class for the module parameter
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef PARAM_H_
#define PARAM_H_

#include <typeinfo>
#include <list>
#include "../assistance/generic_list.h"
#include "param_abstract.h"
#include "param_list.h"

namespace cse_lib {

/// Template class to store parameter for the modules
/**
 * This class defines a data container for the parameter of the modules.  An
 * instance of this class provides an easy access to the data.
 * \tparam is the type of the parameter to store.
 * The Init function can be used to add the instance to a given list of Param
 * instances.
 *
 * \ingroup parameter
 */
template <class T> class Param : public Param_Abstract
{

public:

	/// Create a new instance of the class
	Param()
	{
		config_modified_ = 0;
		configure_by_linked_string_ = false;

		// Automatic string linking in case of bool
		set_bool_str(parameter_);
	}

	virtual ~Param(){}

	/// Getter method with a constant reference as return type
	const T& operator()()
	{
		return parameter_;
	}

	/// Overloaded operator to operate as setter
	void operator() (const T& parameter)
	{
		parameter_ = parameter;
		if(config_modified_ != 0)
			(*config_modified_) = true;
	}

	/// Initialize the parameter class to use in a module
	/**
	 * \param value value to store in the container
	 * \param name name of the container
	 * \param list Parameter list to add the container
	 */
	void Init(const T& value, const string& name, Param_List& list)
	{
		parameter_ = value;
		this->name(name);
		list.Append(this);
		config_modified_ = & (list.config_modified());
	}

	/// Link a certain value to a string
	/**
	 * \param value value to link
	 * \param str string to link with the value
	 */
	void Link_Value_String(T value, string str)
	{
		// Convert string into upper cases
		str = Convert_To_Uppper(str);
		string_list_.Append(str);
		value_list_.Append(value);
		configure_by_linked_string_ = true;
	}

	/// Link the values to the string given in the data structure
	void Link_Value_String(const DS_Enum_String<T>& data_struct)
	{
		for(u_int i = 0; i < data_struct.val1_.length(); i++)
			Link_Value_String(data_struct.val1_[i], data_struct.val2_[i]);
	}

	/// Set the value of the container to the value which is linked to the string
	int Set_Link_Value(string input)
	{
		int index;
		// Convert string into upper cases
		input = Convert_To_Uppper(input);

		index = string_list_.Search(input);
		if(index == -1)
		{
			//Msg(ERROR, "Param", "Instance with name \"" + name_ + "\": "
			//   "No value to given string \"" + input + "\" found, "
			//    "variable keeps his old value: " + this->Str() );
			return -1;
		}

		parameter_ = value_list_.Get_Ref(index);
		if(config_modified_ != 0)
			(*config_modified_) = true;
		return 0;
	}

	/// Return the linked string if available, otherwise return the pure values as string
	string Get_Link_String()
	{
		stringstream temp;

		if(string_list_.length() != 0)
		{
			int index = value_list_.Search(parameter_);
			if(index == -1)
				temp << parameter_;
			else
				temp << string_list_[index];
		}
		else
			temp << parameter_;

		return temp.str();
	}

	/// Set the value(s) of the instance with the data in the given data structure
	void Set_By_Data_Struct(Config_Data_Structure &data)
	{
		string tmp_str;
		T tmp_val;
		int ret;

		// Set value with plain value in data structure or with the linked string
		if(configure_by_linked_string_)
		{
			ret = data.Read_Value(name_, tmp_str, false);

			if(ret != -1)
				if(Set_Link_Value(tmp_str) == -1)
				{
					stringstream s_tmp;
					s_tmp << string_list_;
					throw invalid_argument("\nException[Param]: Configuration of Parameter "
					                       + Quote_Str(name_) + " failed. "
					                       "Given string in configuration " + Quote_Str(tmp_str)
					                       + " does not match to any valid value!\n"
					                       "Valid value are: " + s_tmp.str() + ".");
				}
		}
		else
		{
			try
			{
				ret = Read_Value(data, name_, tmp_val);
			}
			catch(runtime_error&)
			{
				Msg(ERROR, "Param", "Given type \"" + string(typeid(this->parameter_).name()) +
				    "\" of parameter \"" + name_ + "\" cannot be read directly "
				    "out of the configuration structure, "
				    "because no converter function exists. \n"
				    "Solutions: \n "
				    "1) If it is an enum type, please link the possible values with corresponding "
				    "strings in the module's parameter class to solve the problem "
				    "(function:\"Link_Value_String\"\n."
				    "2) Define and implement a \"Read_Value\" function for the given type in the "
				    "\"Config_Data_Strcut\" class!");

				throw;
			}
			catch(invalid_argument& e)
			{
				throw invalid_argument("\nException[Param]: Configuration of Parameter "
				                       + Quote_Str(name_) + " failed." + e.what());
			}

			// In case of a tag was found, write value otherwise the value keeps unchanged.
			if(ret != -1)
				(*this)(tmp_val);
		}

		// Warning if tag was not found and value remains unchanged
		if(ret == -1)
		{
			Msg(WARNING, "Param", "Tag \"" + name_ + "\" not found in current subtree!"
			    " Using default or previous value " + this->Get_Link_String() + ".");

			// Insert value into the given data structure to ensure consistency.
			data.Write_Raw_Value(name_, this->Get_Link_String());
		}
	}

	/// Return the pure stored value as string
	string Str()
	{
		stringstream temp;
		temp << parameter_;
		return temp.str();
	}

	/// Return the status flag, whether the instance can be configured by the linked strings
	bool configure_by_linked_string()
	{
		return configure_by_linked_string_;
	}

	/// Set the status flag, whether the instance can be configured by the linked strings
	void configure_by_linked_string(bool value)
	{
		configure_by_linked_string_ = value;
	}

private:

	virtual void print_parameter(ostream& os)
	{
		// Print the linked string instead of the value if available
		if(string_list_.length() != 0)
		{
			int index = value_list_.Search(parameter_);
			if(index == -1)
				os << name_ << ": "<< parameter_ << "\n";
			else
				os << name_ << ": "<< string_list_[index] << "\n";
		}
		else
			os << name_ << ": "<< parameter_ << "\n";
	}

	// Function for all type, where no special function exists
	template<class F> int Read_Value(Config_Data_Structure &data, std::string name, F& res)
	{
		return data.Read_Value(name, res, false);
	}

	// Function to handle special raw_string data type
	int Read_Value(Config_Data_Structure &data, std::string name, raw_string& res)
	{
		return data.Read_Raw_Value(name, res, false);
	}

	// Special functions and overloading to implement automatic string linking for some types
	void set_bool_str(bool input)
	{
		Link_Value_String(true,"true");
		Link_Value_String(false, "false");
	}

	template<typename TT> void set_bool_str(const TT& input) {}


	// Variable to store the data
	T parameter_;

	// Pointer to the global modified flag
	bool* config_modified_;

	// Lists to assign a string to a certain value
	Generic_List<string, vector> string_list_;
	Generic_List<T, deque> value_list_;

	// flag to store if configuration is done with the linked string or not
	bool configure_by_linked_string_;
};

/// Template class to store the calculated parameters of the parameter units.
/**
 * This class defines a data container for the calculated parameter.  An
 * instance of this class provides an easy access to the data.
 * \tparam is the type of the parameter to store.
 * The Register function can be used to add the instance to a given list of Param_Out
 * instances.
 *
 * \ingroup parameter
 */
template <class T> class Param_Out : public Param_Abstract
{

public:

	/// Create a new instance of the class
	Param_Out()
	{
		once_written_ = false;
		// Automatic string linking in case of bool
		set_bool_str(parameter_);
	}

	/// Create a new instance and set the value to init_value
	explicit Param_Out(const T init_value)
	{
		parameter_ = init_value;
		once_written_ = true;
		// Automatic string linking in case of bool
		set_bool_str(parameter_);
	}

	virtual ~Param_Out(){}

	/// Overloaded operator to operate as getter
	/*T operator()() const
	{
		return parameter_;
	}*/

	T& operator()()
	{
		return parameter_;
	}

	/// Overloaded operator to operate as setter
	void operator() (T parameter)
	{
		parameter_ = parameter;
		once_written_ = true;
	}

	/// Register the parameter instance to use in a module
	/**
	 * \param name Name of the parameter output port
	 * \param list List to register the port
	 */
	void Register(string name, Param_List &list)
	{
		this->name(name);
		list.Append(this);
	}

	/// Link a value with a string
	/**
	 * \param value value to link
	 * \param str string to link with the value
	 */
	void Link_Value_String(T value, string str)
	{
		// Convert string into upper cases
		str = Convert_To_Uppper(str);
		string_list_.Append(str);
		value_list_.Append(value);
	}

	/// Link the values to the string given in the data structure
	void Link_Value_String(DS_Enum_String<T> data_struct)
	{
		for(u_int i = 0; i < data_struct.val1_.length(); i++)
			Link_Value_String(data_struct.val1_[i], data_struct.val2_[i]);
	}

	/// Set the value of the container to the value which is linked to the string
	int Set_Link_Value(string input)
	{
		int index;
		// Convert string into upper cases
		input = Convert_To_Uppper(input);

		index = string_list_.Search(input);
		if(index == -1)
		{
			Msg(ERROR, "Param_Out", "Instance with name \"" + name_ + "\": "
			    "No value to given string \"" + input + "\" found, "
			    "variable keeps his old value: " + this->Str() );
			return -1;
		}

		parameter_ = value_list_[index];
		once_written_ = true;

		return 0;
	}

	/// Return the linked string if available, otherwise return the pure values as string
	string Get_Link_String()
	{
		stringstream temp;

		if(string_list_.length() != 0)
		{
			int index = value_list_.Search(parameter_);
			if(index == -1)
				temp << parameter_;
			else
				temp << string_list_[index];
		}
		else
			temp << parameter_;

		return temp.str();
	}

	/// Return the pure stored value as string
	string Str()
	{
		stringstream temp;
		temp << parameter_;
		return temp.str();
	}

private:

	virtual void print_parameter(ostream& os)
	{
		// Print the linked string instead of the value if available.
		if(string_list_.length() != 0)
		{
			int index = value_list_.Search(parameter_);
			if(index == -1)
				os << name_ << ": "<< parameter_ << "\n";
			else
				os << name_ << ": "<< string_list_[index] << "\n";
		}
		else
			os << name_ << ": "<< parameter_ << "\n";
	}

	// Special functions and overloading to implement automatic string linking for some types
	void set_bool_str(bool input)
	{
		Link_Value_String(true,"true");
		Link_Value_String(false, "false");
	}

	template<typename TT> void set_bool_str(const TT& input) {}

	/// Variable to store the data
	T parameter_;

	/// Flag to indicate, when the parameter has been written (and is valid).
	bool once_written_;

	// Lists to assign a string to a certain value
	Generic_List<string, vector> string_list_;
	Generic_List<T, deque> value_list_;
};

}
#endif


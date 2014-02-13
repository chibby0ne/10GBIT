//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Abstract base class for the status out port
/// \author Timo Lehnigk-Emden
/// \date   2011/01/14
//

#ifndef STATUS_OUT_ABSTRACT_H_
#define STATUS_OUT_ABSTRACT_H_

#include <vector>
#include "../../assistance/generic_list.h"

namespace cse_lib {


/// Data structure to store the status out value
/**
 * This data structure hold a output value of the status out class for further
 * processing.  Each value as one position, consisting of any number of
 * addresses.
 */
class Status_Out_Data_Struct
{
public:

	/// Value as string
	string value;

	/// List contains one address value for each dimension
	Generic_List<string, vector> address;

	/// List contains a name for each dimension
	Generic_List<string, vector> dim_name;

	/// Delete all stored values of the data container
	void Clear()
	{
		value.clear();
		address.Clear();
		dim_name.Clear();
	}
};

// Forward declartion of the list to use it as parameter of a member function
class Status_Out_List;

/// Abstract class for the status port of the modules
/**
 * \ingroup base_class
 */
class Status_Out_Abstract
{
public:
	Status_Out_Abstract()
	{
		write_into_results_ = true;
	}

	virtual ~Status_Out_Abstract() {}

	/// Abstract function, which is overwritten in the derived classes
	virtual bool Get_Value_Str(Status_Out_Data_Struct& data) = 0;

	/// Return the name of the status out port
	string name()
	{
		return name_;
	}

	/// Set the name of the status out port
	void name(string name)
	{
		name_ = name;
	}

	/// Getter for write_into_results member variable
	/**
	 * If true (=default) the result stored in the ports can be read by the Get_Value-Str() function
	 */
	bool write_into_results()
	{
		return write_into_results_;
	}

	/// Set for write_into_results
	void write_into_results(bool value)
	{
		write_into_results_ = value;
	}

	/// Register the port into a list and set a name
	/**
	* \param name name of the port
	* \param list list to add the port
	* \param write_into_results  = true (=default) the results stored in this port is written
	* into a data structure (e.g., XML result file) if the Managa_Module_Results class is used
	*/
	void Register(string name, Status_Out_List &list, bool write_into_results = true);

protected:

	/// Name of the container
	string name_;

	/// Should this port automatically read and written to result data
	bool write_into_results_;
};

}
#endif

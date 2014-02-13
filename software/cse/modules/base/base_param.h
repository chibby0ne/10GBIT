//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Abstract interface class for parameter of the modules
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef BASE_PARAM_H_
#define BASE_PARAM_H_


#include "cse_namespace.h"
#include "../assistance/config_data_struct.h"
#include "param_list.h"
#include "param.h"

namespace cse_lib {

/// Base of the parameter classes for the modules
/**
 * \ingroup base_class
 */
class Base_Parameter
{

public:
	virtual ~Base_Parameter() {}

	/// Instance name parameter
	Param<string> instance_name;

	/// Read parameter from a data tree structure.
	/**
	 * \param config_tree Class from which to read the configuration parameter.
	 */
	virtual void Config_From_Data_Struct(Config_Data_Structure& config_tree)
	{
		// Iterate over all parameters in the list and call the configure function
		Param_List::iterator& iter = param_list_.iter;

		for(iter = param_list_.begin(); iter != param_list_.end(); iter++)
			(*iter)->Set_By_Data_Struct(config_tree);
	};

	/// Return a reference to the parameter list object
	Param_List& param_list()
	{
		return param_list_;
	}

	/// Return a reference to the output parameter list object
	Param_List& param_out_list()
	{
		return param_out_list_;
	}

protected:

	/// List to store all parameter objects
	Param_List param_list_;

	/// List to store all output parameter object (param units only)
	Param_List param_out_list_;

	/// Overload streaming operator, streams all objects which are in the parameter list param_list_.
	inline friend ostream & operator<<(ostream & os, Base_Parameter & data)
	{
		os << "\nList of all parameters of module "<< data.instance_name() <<":\n" << data.param_list_ << endl;
		return os;
	}
};
}
#endif

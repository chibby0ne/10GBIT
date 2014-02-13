//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Parameter template class for the module parameter definition
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef PARAM_ABSTRACT_H_
#define PARAM_ABSTRACT_H_

#include <iostream>
#include <stdexcept>

namespace cse_lib {

/// Abstract class of class Param
/**
 * Only to provide a common type to store parameter classes of different types
 * in a list.  In addition a common streaming operator is defined.
 *
 * \ingroup base_class
 */
class Param_Abstract
{

public:

	Param_Abstract(){}
	virtual ~Param_Abstract(){}

	/// Streams the stored values in the container
	inline friend ostream & operator<<(ostream & os, Param_Abstract & data)
	{
		data.print_parameter(os);
		return os;
	}

	/// Return the name of the parameter class
	string name()
	{
		return name_;
	}

	/// Set the name of the parameter class
	void name(string name)
	{
		name_ = name;
	}

	/// Virtual function to get access to the values of the derived classes (interface function)
	virtual string Str() = 0;

	/// Virtual function to get access to the values of the derived classes (interface function)
	virtual string Get_Link_String() = 0;

	/// Virtual function to get access to the values of the derived classes (interface function)
	virtual void Set_By_Data_Struct(Config_Data_Structure &data)
	{
		Msg(ERROR, "Param_Abstract", "Instance with name \"" + name() + "\": "
		    "The function Set_By_Data_Struct has no functionality for this module! ");
		throw runtime_error("Programming error! Please report this bug to the developer team");
	}

protected:

	/// Name of the container
	string name_;

private:

	virtual void print_parameter(ostream &os) = 0;
};
}
#endif

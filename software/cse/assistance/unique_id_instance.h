//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Create an instance of a module with the given unique ID
/// \author Timo Lehnigk-Emden
/// \date   2010/12/14
//
#ifndef UNIQUE_ID_INSTANCE_H_
#define UNIQUE_ID_INSTANCE_H_

#include <string>
#include <exception>
#include "../modules/base/base_param.h"
#include "generic_list.h"
#include "helper_functions.h"

/// Abstract class for the template class Unit_Instance_Manager
/**
 * \ingroup base_class
 */
class Unit_Instance_Manager_Abstract
{
public:
	virtual ~Unit_Instance_Manager_Abstract(){};

	/// Abstract function interface
	virtual cse_lib::Base_Parameter* Create_Instance(std::string unique_id)=0;
};

/// Class to create a instance of a module or unit dynamically by its unique id
/**
 * This template class provides for each module type an function to create an instance of this
 * type dynamically. \tparam MODULE_TYPE is the type of the module or unit
 */
template<class MODULE_TYPE> class Unit_Instance_Manager : public Unit_Instance_Manager_Abstract
{
public:
	/// Create a new instance and add the instance to the given list
	Unit_Instance_Manager(Generic_Pointer_List<Unit_Instance_Manager_Abstract,std::vector > list)
	{
		list.Append(this);
	}

	Unit_Instance_Manager(){}
	virtual ~Unit_Instance_Manager(){}

	/// Create an instance of type MODULE_TYPE
	/**
	 * The instance is created only if the unique id matches of this type
	 * matches the given unique_id
	 * \param unique_id Unique id of the instance to create
	 */
	cse_lib::Base_Parameter* Create_Instance(std::string unique_id)
	{
		cse_lib::Base_Parameter* ptr = 0;
		if (unique_id == MODULE_TYPE::Unique_ID())
		{
			try
			{
				ptr = new MODULE_TYPE;
			}
			catch (std::bad_alloc&)
			{
				hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Unit_Instance_Manager",
				    "Unit or module with the unique ID \"" + unique_id + "\" "
				    "cannot be instantiated, because of a memory allocation failure!");

				throw;
			}
		}
		return ptr;
	}
};
#endif

//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  
/// \author Timo Lehnigk-Emden
/// \date   2010/12/15
//

#ifndef MANAGE_DYN_INSTANCES_H_
#define MANAGE_DYN_INSTANCES_H_

#include "unique_id_instance.h"
#include "../modules/base/base_param.h"

/// Mangage the dynamic instances of modules or parameter units
/**
 * \ingroup base_class
 */
class Manage_Dynamic_Instances
{
public:
	Manage_Dynamic_Instances()
	{
		// Register default instance only once.
		if(inst_cnt_ == 0)
			Make_Modules_Visible();

		inst_cnt_++;
	}

	virtual ~Manage_Dynamic_Instances()
	{
		Delete_All_Modules();

		// Delete all objects from static list if this is the last instance
		if(inst_cnt_ == 1)
			Unregister_All_Create_Modules();

		inst_cnt_--;
	}

	/// Create an instance of a module with a given unique name and returns its pointer
	cse_lib::Base_Parameter* Create_Instance_By_Name(std::string unique_id)
	{
		cse_lib::Base_Parameter* module_ptr = 0;
		cse_lib::Base_Parameter* module_ptr_tmp = 0;

		// Iterator over the list
		Generic_Pointer_List<Unit_Instance_Manager_Abstract,std::vector >::iterator iter;

		// Go thru all objects
		for(iter = unit_inst_list_.begin(); iter != unit_inst_list_.end(); iter++)
		{
			module_ptr_tmp = (*iter)->Create_Instance(unique_id);

			if(module_ptr_tmp != 0)
			{
				iter++;
				break;
			}
		}
		module_ptr = module_ptr_tmp;

		// Report if no module was found
		if(module_ptr_tmp == 0)
			hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Manage_Dynamic_Instances",
			                    "No object of type \"Unit_Instance_Manger\" to create a module or unit"
		                        "with the given unique id \"" + unique_id + "\" was found in the list!");

		// Search for double entries
		for(; iter != unit_inst_list_.end(); iter++)
		{
			module_ptr_tmp = (*iter)->Create_Instance(unique_id);
			if(module_ptr_tmp != 0)
				hlp_fct::logging::Msg(hlp_fct::logging::WARNING, "Manage_Dynamic_Instances",
				                    "The object of type \"Unit_Instance_Manger\" to create a module or unit"
				                    "with the given unique id \"" + unique_id + "\" is twice in the list."
				                    "Please report this message to the developer!");
		}

		return module_ptr;
	}

	/// Add a new instance with the given unique ID to the list
	cse_lib::Base_Parameter* Add_Instance_By_Name(std::string unique_id, bool status_out = false)
	{
		cse_lib::Base_Parameter* ptr = 0;
		ptr = Create_Instance_By_Name(unique_id);
		if(ptr != 0)
		{
			unit_list_.Append(ptr);
			if(status_out)
				hlp_fct::logging::Msg(hlp_fct::logging::INFO, "Manage_Dynamic_Instances",
				                    "New instance of object \"" + unique_id + "\""
				                    "created successfully!");
		}
		return ptr;
	}

	/// Search a module / unit by its instance name
	cse_lib::Base_Parameter* Search_by_Instance_Name(std::string name, bool error_out = true)
	{
		cse_lib::Base_Parameter* ptr = 0;

		for(unit_iter_ = unit_list_.begin(); unit_iter_ != unit_list().end(); unit_iter_++)
		{
			if((*unit_iter_)->instance_name() == name)
			{
				ptr = (*unit_iter_);
				break;
			}
		}

		// Report if no module was found
		if(ptr == 0 && error_out)
			hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Manage_Dynamic_Instances", "No module or unit"
			                    "with the given instance name \"" + name + "\" was found in the list!");

		return ptr;
	}

	/// Delete all instantiated modules in the list
	void Delete_All_Modules()
	{
		// Delete all objects
		for(unit_iter_ = unit_list_.begin(); unit_iter_ != unit_list_.end(); unit_iter_++)
			delete (*unit_iter_);

		// Delete all list entries
		unit_list_.Clear();
	}

	/// Returns a reference to the list which holds all pointer of the instantiated modules
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>& unit_list()
	{
		return unit_list_;
	}

	/// Returns a reference to the iterator for the pointer list unit_list
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>::iterator& unit_iter()
	{
		return unit_iter_;
	}

	// Helper function to have an easy syntax to register a module
	template<class T> static void Register_Module()
	{
		unit_inst_list_.Append(new Unit_Instance_Manager<T>);
	}
private:
	// List to store the instance creation objects for each module
	static Generic_Pointer_List<Unit_Instance_Manager_Abstract, std::vector> unit_inst_list_;

	// Count the number created object instances to handle the static list unit_inst_list_
	static unsigned int inst_cnt_;

	/// List to store param units or modules
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector> unit_list_;

	// Iterator for the unit_list_
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>::iterator unit_iter_;

	// Delete all objects, which can create a module from the statit list unit_inst_list_
	void Unregister_All_Create_Modules()
	{
		// Delete all objects with can dynamic instantiate modules
		for(unsigned int i = 0; i < unit_inst_list_.length();i++)
			delete unit_inst_list_.Get(i);
	}

	// Delete all list objects
	/*void Unregister_All_Modules()
	{
		// Delete all dynamic created modules / units
		for(unit_iter_ = unit_list_.begin(); unit_iter_ != unit_list_.end(); unit_iter_++)
			delete (*unit_iter_);

	}*/

	// Make all modules visible for an automatic instantiation
	void Make_Modules_Visible();
};

#endif

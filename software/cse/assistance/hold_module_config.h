//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Holds the configuration of the given modules and allows the automatic configuration of the modules.
/// \author Timo Lehnigk-Emden
/// \date   2011/11/25
//

#ifndef HOLD_MODULE_CONFIG_H_
#define HOLD_MODULE_CONFIG_H_

#include "config_data_struct_pugixml.h"
#include "generic_pointer_list.h"
#include "../modules/base/base_param.h"

/// This class hold one configuration set for the modules.
/**
 * You can copy this class, without copy the configuration sets, this allows for an memory efficient
 * storage of the currently used configuration. This class does not manage the allocation and deletion
 * of the objects behind the stored pointer. The memory management must be done outside. This implementation
 * allows to copy the class without unneeded copy operations and memory efficient storage.
 */
class Hold_Module_Config
{
public:

	Hold_Module_Config() {}

	virtual ~Hold_Module_Config() {}

	/// Configure a module with the section in the XML file.
	/**
	 * \param module Module to configure with the stored parameters.
	 */
	Config_Data_Struct_Pugixml* Configure_Module(cse_lib::Base_Parameter &module) const;

	/// Type for the module configuration list.
	typedef Generic_Pointer_List<Config_Data_Struct_Pugixml, ::std::vector> Xml_Module_List;

	/// Type for the list iterator.
	typedef Xml_Module_List::iterator Xml_Module_List_Iter;

	/// Returns a reference to the module list.
	Xml_Module_List& xml_module_list()
	{
		return xml_module_list_;
	}

	/// Compare operator
	bool operator== (const Hold_Module_Config& rhs )
	{
		if(xml_module_list_ == xml_module_list_)
			return true;
		else
			return false;
	}

	/// Clear and delete all objects in the list
	/**
	 * This function execute a deletion to all pointers which are stored in the list and set the
	 *  list size to zero.
	 */
	void Destroy_List()
	{
		Xml_Module_List_Iter iter;

		for(iter = xml_module_list_.begin(); iter != xml_module_list_.end(); iter++)
			delete (*iter);

		xml_module_list_.Clear();
	}

	/// Streaming operator to print out the current module configuration set
	inline friend std::ostream & operator<<(std::ostream & os, Hold_Module_Config& cfg)
	{
		Xml_Module_List_Iter iter;

		for(iter = cfg.xml_module_list_.begin(); iter != cfg.xml_module_list_.end(); iter++)
			os << *(*iter);

		return os;
	}

private:

	// Find a xml module section from the list by the instance name and returns a pointer the result
	Config_Data_Struct_Pugixml* Find_XML_Object_By_Name(const std::string& name) const;

	// List and iterator to store XML subtrees for each module
	mutable Xml_Module_List xml_module_list_;

	// Fixed class name as string
	static ::std::string instance_name(){ return "Hold_Module_Config";}
};


#endif

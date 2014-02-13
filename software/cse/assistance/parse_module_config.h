//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Implementation of the parse configuration of the modules from XML file
/// \author Timo Lehnigk-Emden
/// \date   2010/11/09
//

#ifndef PARSE_MODULE_CONFIG_H_
#define PARSE_MODULE_CONFIG_H_

#include <string>
#include <vector>
#include <iostream>
#include "config_data_struct_pugixml.h"
#include "generic_pointer_list.h"
#include "manage_dyn_instances.h"
#include "hold_module_config.h"

/// Parse the configuration XML file to configure the module instances.
/** 
 * This class parsing a XML configuration file and search for the "module" sections.
 * It provides a function the configure a module with the parameter given in the XML file.
 * The assignment of the XML section is done over the instance name!
 * In addition the module automatic instaziate the given parameter unit modules and propagate
 * their ports into to global variable section. After that the global variables were propagated into
 * the module sections.
 *
 * \ingroup parameter
 */
class Parse_Module_Config : public Config_Data_Struct_Pugixml
{
public:

	/// Create a new instance of the class, use the Reload function to read in XML data.
	Parse_Module_Config(int max_number_print_into_tag = -1) : MAX_NUMBER_PRINT_INTO_TAG(max_number_print_into_tag) {}

	/// Parse a XML configuration file.
	/**
	 * \param file_name file name and path of the XML file.
	 * \param max_number_print_into_tag the number of elements to write into the result file, if
	 *  the parameter (xml tag) is a list of values.
	 *  This feature is useful to keep the result files small and readable.
	 */
	explicit Parse_Module_Config(::std::string file_name, int max_number_print_into_tag = -1);

	/// Parse a XML stream
	/**
	 * \param stream C++ stream which contains the XML tree.
	 * \param max_number_print_into_tag the number of elements to write into the result file, if
	 *  the parameter (xml tag) is a list of values.
	 *  This feature is useful to keep the result files small and readable.
	 */
	explicit Parse_Module_Config(::std::iostream& stream, int max_number_print_into_tag = -1);

	virtual ~Parse_Module_Config();

	/// Configure a module with the section in the XML file.
	/**
	 * \param module Module to configure with the parameter in the XML file.
	 */
	void Configure_Module(cse_lib::Base_Parameter &module) const;


	/// Reload the XML data from a c++ stream.
	/**
	 * \param stream c++ stream which contains the new XML tree.
	 */
	void Reload(::std::iostream& stream);

	/// Reload the XML data form a XML document (copy the given document tree).
	/**
	 * \param doc xml document
	 */
	void Reload(const pugi::xml_document& doc);

	/// Returns a reference of the current XML tree.
	pugi::xml_node Get_Current_Tree()
	{
		return doc_.document_element();
	}

	/// Return the list which contains all instantiated parameter units.
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>& Get_Param_Unit_List()
	{
		return dyn_inst_.unit_list();
	}

	/// Returns the object which hold the current configuration for the modules.
	const Hold_Module_Config& Get_Current_Module_Config()
	{
		return current_module_config_;
	}

	/// Print the configuration XML tree of the given module as template for the configuration file.
	static void Print_Module_Config_XML_Tree(cse_lib::Base_Parameter& module);

	/// String contains the name of the attribute to mark the nodes which can be shortened in the XML file.
	static const std::string SHORTENED_LIST_ATTRIBUTE;

	/// String contains the prefix for all internal attributes.
	static const std::string INTERNAL_ATTRIBUTE_PREFIX;

private:

	// Find a global variable in the XML tree
	pugi::xml_node Find_Global_Variable(std::string name, std::string& value);

	// Replace the global variables references in the XML tree
	void Replace_Global_Variable(pugi::xpath_node);

	// Read the module sesctions
	void Read_Module_Sections();

	// Read the parameter section and process them
	void Read_Param_Units();

	// Delete all objects in the xml_module_list_ and the list itself
	void Clear_List();

	bool Find_Param_Unit_Variable(std::string unit_name, std::string var_name,
	                              std::string& param_value);

	// Create an instance of the object with the given section and configure it
	void Create_Instance_and_Configure(pugi::xpath_node node);

	void Process_Param_Units();

	// Execute the Run function of the given object
	void Exec_Param_Unit(cse_lib::Base_Parameter* param_unit);

	void Replace_Processed_Param_Ports(pugi::xpath_node node, std::string instance);

	// Insert the updated missing parameter of a module after the configuration
	void Insert_Missing_Module_Parameter(Config_Data_Struct_Pugixml& module) const;

	// Insert all output parameter and values from the given parameter unit into the current XML tree.
	void Insert_Param_Unit_Output_Values(cse_lib::Base_Parameter& unit);

	// Check if a given parameter port exist in a module or unit
	cse_lib::Param_Abstract& Check_Param_Port(cse_lib::Base_Parameter& obj, ::std::string port_name, bool status_out = true);

	// Check weather the node has a CSE internal attribute
	bool Detect_Internal_Attribute(pugi::xml_node node);

	// Copy internal attribute from src to dst
	void Copy_Internal_Attributes(pugi::xml_node src, pugi::xml_node dst);

	static ::std::string instance_name(){ return "Parse_Module_Config";}

	// Instance to create dynamically an instance of a module or unit
	Manage_Dynamic_Instances dyn_inst_;

	// Internal variable to store the number of elements of the list tags to print into the result file.
	int MAX_NUMBER_PRINT_INTO_TAG;

	// Holds to current configuration of the modules (pointer list)
	Hold_Module_Config current_module_config_;
};


#endif

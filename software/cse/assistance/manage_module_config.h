//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  
/// \author Timo Lehnigk-Emden
/// \date   2010/12/20
//

#ifndef MANAGE_MODULE_CONFIG_H_
#define MANAGE_MODULE_CONFIG_H_

#include <string>
#include "parse_module_config.h"
#include "config_data_struct_pugixml.h"
#include "iter_var_struct.h"
#include "../modules/base/base_param.h"

/// Manage the automatic configuration of the modules via a XML file.
/**
 * This class enhance the configuration via the XML file with iteration variables in the XML file.
 * That means several different configuration can be processed with the variables in the global/iter
 * section of the XML file. This class uses the Parse_Module_Config class.
 *
 * \ingroup infrastructure
 */
class Manage_Module_Config : public Config_Data_Struct_Pugixml
{
public:

	/// Create a new instance for automatic configuration via a XML file
	/**
	 * \param filename XML file to use for read the configuration.
	 * \param max_number_print_into_tag the number of elements to write into the result file, if
	 *  the parameter (xml tag) is a list of values.
	 *  This feature is useful to keep the result files small and readable.
	 */
	explicit Manage_Module_Config(std::string filename, int max_number_print_into_tag = -1);

	/// Create a new instance for automatic configuration via XML data stored in a stream
	/**
	 * \param stream C++ stream contains the XML data to use for read the configuration.
	 * \param max_number_print_into_tag the number of elements to write into the result file, if
	 *  the parameter (xml tag) is a list of values.
	 *  This feature is useful to keep the result files small and readable.
	 */
	explicit Manage_Module_Config(std::iostream& stream, int max_number_print_into_tag = -1);


	virtual ~Manage_Module_Config();


	/// Configure the given module with the data from the XML file
	/** The function searches in the XML file a module section where the instance_name tag matches
	 * the instance name of the given module. The function throw an exception if no matching
	 * XML section was found
	 */
	void Configure_Module(cse_lib::Base_Parameter& module);

	/// Update configuration instance with the next iteration point
	/**
	 * \return -1 if all variables have been iterated or no iteration variable
	 * is present, otherwise 0
	 */
	int Update_To_Next_Iter();


	/// Return the list which contains all instantiated parameter units
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>& Get_Param_Unit_List()
	{
		return curr_config_.Get_Param_Unit_List();
	}

	/// Type of list which holds the iteration variable objects
	typedef Generic_List<Iteration_Variable_Structure, std::vector> ITER_VAR_LIST_TYPE;

	/// Return the list which contains all iteration variable object
	ITER_VAR_LIST_TYPE& iter_var_list()
	{
		return data_set_list_.Get_Ref(curr_data_set_index_);
	}

	/// Return a reference of the current configuration object
	Parse_Module_Config& Get_Current_Config()
	{
		return curr_config_;
	}

	/// Returns a reference of the current XML tree
	pugi::xml_node Get_Inital_Tree()
	{
		return initial_doc_.document_element();
	}

	/// Streaming operator to print out the current XML tree
	inline friend std::ostream & operator<<(std::ostream & os, Manage_Module_Config& data)
	{
		os << data.curr_config_;
		return os;
	}

	/// Write the requested parameter by module name and parameter name in the given variable
	/**
	 * Only parameters form param_units which are instantiated automatically can be requested
	 * \param unit_name Name of the unit which contains the requested parameter
	 * \param param_name Name of the parameter
	 * \param value Variable to write the value of the correct type
	 * \return 0 in case of no error
	 *         1 Type miss-match between given variable type and type of requested parameter
	 *         2 Variable cannot be found (wrong module or parameter name)
	 */
	template <class T> int Get_Parameter_By_Name(std::string unit_name,
	                                             std::string param_name,
	                                             T& value) const;

	/// Register a new module or parameter unit
	template<class T> static void Register_Module()
	{
		Manage_Dynamic_Instances::Register_Module<T>();
	}

private:

	// Read iteration variable
	void Read_Iter_Var();

	// Create configuration instances
	void Update_Config_Instances(ITER_VAR_LIST_TYPE& var_list);

	// Replace iteration variable with var_name by value
	void Replace_Iter_Var(std::string var_name, std::string value);

	typedef Generic_List<std::string, std::vector> STRING_LIST_TYPE;

	// Process iteration loops, process and store the variable values and names
	void Process_Iter_Loops(pugi::xml_node& node,
	                        ITER_VAR_LIST_TYPE& var_list,
	                        STRING_LIST_TYPE& name_list);

	// Process iteration data set
	void Process_Iter_Data_Sets();

	// Update the iteration variable of a data set
	int Next_Iter(ITER_VAR_LIST_TYPE& var_list);

	typedef Generic_List<ITER_VAR_LIST_TYPE, std::vector> DATA_SET_LIST_TYPE;

	// Update the iteration variable to the next data set
	int Next_Data_Set();

	//ITER_VAR_LIST_TYPE var_list_;
	ITER_VAR_LIST_TYPE::iterator var_iter_;

	// List to store the data sets lists
	DATA_SET_LIST_TYPE data_set_list_;

	// Object hold the current valid configuration
	mutable Parse_Module_Config curr_config_;

	pugi::xml_document initial_doc_;

	unsigned int curr_data_set_index_;
};

// Template function implementation because of linking issues
template <class T> int Manage_Module_Config::Get_Parameter_By_Name(std::string param_unit,
                                                                   std::string param_name,
                                                                   T& value) const
{
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>::iterator mod_iter;
	Generic_Pointer_List<cse_lib::Base_Parameter, ::std::vector>&
	list = curr_config_.Get_Param_Unit_List();
	bool found_in_out = false;

	cse_lib::Base_Parameter* mod_ptr = 0;

	// Find module
	for(mod_iter = list.begin(); mod_iter != list.end(); mod_iter++)
	{
		if((*mod_iter)->instance_name() == param_unit)
		{
			mod_ptr = (*mod_iter);
			break;
		}
		std::cout << (*mod_iter)->instance_name() << std::endl;
	}

	if(mod_ptr == 0)
	{
		hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Manage_Module_Config",
		                      "No module or unit with name \"" + param_unit + "\" has been found!");

		hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Manage_Module_Config","Available Units are: ", false);
		for(mod_iter = list.begin(); mod_iter != list.end(); mod_iter++)
			Msg(hlp_fct::logging::ERROR, "Manage_Module_Config", (*mod_iter)->instance_name(), false);
		return 2;
	}

	// Find parameter in input parameter
	cse_lib::Param_Abstract* param_ptr = mod_ptr->param_list().Search(param_name);

	// If not found check output parameter
	if(param_ptr == 0)
	{
		param_ptr = mod_ptr->param_out_list().Search(param_name);
		found_in_out = true;
	}

	// Check weather parameter exists
	if(param_ptr == 0)
	{
		hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Manage_Module_Config",
		                      "Requested parameter \"" + param_name + "\" "
		                      "of module \"" + param_unit + "\" does not exits!");
		return 2;
	}

	// Check for right type
	cse_lib::Param<T>* value_ptr = dynamic_cast<cse_lib::Param<T>* > (param_ptr);
	cse_lib::Param_Out<T>* value_out_ptr = dynamic_cast<cse_lib::Param_Out<T>* > (param_ptr);

	if(value_ptr == 0 && value_out_ptr == 0)
	{
		hlp_fct::logging::Msg(hlp_fct::logging::ERROR, "Manage_Module_Config",
		                      "Requested parameter \"" + param_name + "\" of module \""
		                      + param_unit + "\" has not the same type than the given target variable!");
		return -1;
	}

	// Choose right pointer
	if(found_in_out)
		value = (*value_out_ptr)();
	else
		value = (*value_ptr)();

	return 0;
}

#endif

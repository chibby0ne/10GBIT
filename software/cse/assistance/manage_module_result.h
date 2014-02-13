//
//  Copyright © 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems
//
/// \file
/// \brief  Manage the result file for the modules / chains
/// \author Timo Lehnigk-Emden
/// \date   2011/01/13
//

#ifndef MANAGE_MODULE_RESULT_H_
#define MANAGE_MODULE_RESULT_H_

#include <string>
#include "../modules/base/base_iface.h"
#include "result_data_struct_pugixml.h"
#include "parse_module_config.h"
#include "iter_var_struct.h"
#include "manage_module_config.h"


/// Manage the XML result file
/**
 * \ingroup infrastructure
 */
class Manage_Module_Result
{
public:

	/// Create a new XML result file
	explicit Manage_Module_Result(std::string filename = "last_result.xml",
	                              std::string executable_name = "unknown");

	virtual ~Manage_Module_Result();

	/// Write the current XML result tree into the file
	void Write_Current_State();

	/// Insert the result from a module into the current working tree
	void Insert_Results_From_Module(cse_lib::Base_Interface& module);

	/// Create a new iteration value XML tree to store the modules results
	/**
	 * Create a new subtree which contains the current iteration variable states
	 * and the current configuration for all modules!
	 * \param config object contains the configuration of the modules
	 * and the iteration variable objects
	 */
	void Create_Iteration_Value_Result_Point(Manage_Module_Config& config);

	/// Set the free comment tag in the meta information section.
	void Set_Free_Comment(const std::string& text);

	/// Append a text to the free comment tag in the meta information section.
	void Append_Free_Comment(const std::string& text);

private:

	typedef Generic_List<Iteration_Variable_Structure, std::vector> ITER_VAR_LIST_TYPE;

	/// Insert next iteration variable
	void Insert_Iteration_Point(ITER_VAR_LIST_TYPE list);

	/// Insert configuration into the current subtree
	void Insert_Configuration(Parse_Module_Config& curr_config);

	/// Insert initial configuration into the file
	void Insert_Initial_Configuration(Manage_Module_Config& config);

	/// Insert a set of meta information on the simulation in the result file
	void Insert_Meta_Information();

	// XML file object
	Result_Data_Struct_Pugixml xml_file_;

	struct
	{
		/// Meta information executable name
		std::string executable_name;

		/// Meta information simulation date
		std::string execution_date;

		/// Free text, comment
		std::string free_comment;

	} meta_information_;
};

#endif

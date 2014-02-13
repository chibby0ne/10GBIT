//
//  Copyright © 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems
//
/// \file
/// \brief  Result data struct with the pugi-xml parser to write XML files
/// \author Timo Lehnigk-Emden
/// \date   2011/01/13
//

#ifndef RESULT_DATA_STRUCT_PUGIXML_H_
#define RESULT_DATA_STRUCT_PUGIXML_H_

#include "result_data_struct.h"
#include "../third_party/pugixml/pugixml.hpp"

/// Creates a new XML file and writes tags and attributes with the pugixml library
/**
 * \ingroup infrastructure
 */
class Result_Data_Struct_Pugixml : public Result_Data_Struct
{
public:
	/// Create a new empty XML file
	/**
	 * \param filename name of the XML file
	 * \param top_tag unique top level tag
	 */
	explicit Result_Data_Struct_Pugixml(std::string filename, std::string top_tag);

	virtual ~Result_Data_Struct_Pugixml();

	/// Set the current tree pointer to the root of the tree
	void Goto_Root();

	/// Set the current tree pointer to tree with the name param_name
	int Goto_Tree(::std::string tree_name, bool error_out=true);

	/// Set the current tree pointer to a subtree with a name
	int Goto_Sub_Tree(::std::string tree_name, bool error_out=true);

	/// Go up in the tree
	/**
	* \param levels number of level to go up (default = 1)
	*/
	void Goto_Up(unsigned int levels = 1);

	/// Create a new tag in the XML file
	void Create_Tag(::std::string tag_name);

	/// Write new data
	void Create_Value(::std::string branch_name, ::std::string value);

	/// Write attribute to the current value
	void Create_Attribute(::std::string attrib_name, ::std::string value);

	/// Delete the child tag with the given name
	void Delete_Tag(const ::std::string& tag_name);

	/// Save the current state to the file
	void Save_File();

	/// Write a given pugixml subtree at the current tree position
	void Insert_Subtree(pugi::xml_node new_subtre);

	/// Store the current position of the pointer in the XML tree
	void Save_Tree_Position()
	{
		save_pos_ = current_node_;
	}

	/// Restore the previously saved tree position
	void Restore_Tree_Position()
	{
		current_node_ = save_pos_;
	}

	/// Short the values of all the nodes in the current sub tree, which are marked with a special attribute
	int Short_Marked_Nodes();

private:
	/// Variable to hold the xml document
	pugi::xml_document doc_;

	/// Current node for processing
	pugi::xml_node current_node_;

	/// Filename and path of the XML file
	std::string filename_;


	/// Check for multiple node with the same name
	int Check_Multiple_Tags(pugi::xpath_node_set node_set, ::std::string tag_name,
	                        bool status_out=true);

	/// Search for a given node in the subtree of the given node and returns a list of nodes
	pugi::xpath_node_set Search_Node(pugi::xml_node node, ::std::string tag_name,
	                                 bool status_out=true);

	/// Short the given comma separated list after the number-th comma
	void Short_String(std::string& var_value, int number);

	// Saved position
	pugi::xml_node save_pos_;
};

#endif

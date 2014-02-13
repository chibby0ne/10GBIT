//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the xml class which uses the pugixml library
/// \author Timo Lehnigk-Emden
/// \date   2010/11/11
//

#ifndef CONFIG_DATA_STRUCT_PUGIXML_H_
#define CONFIG_DATA_STRUCT_PUGIXML_H_

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "config_data_struct.h"
#include "../third_party/pugixml/pugixml.hpp"
#include "generic_list.h"

/// XML Parser interface for the pugixml parser
/** This class is one implementation of the abtract Config_Data_Struct interface class
 * It provides a simple interface for the PugiXML Parser
 * \ingroup assistance
 * \ingroup parameter
 */
class Config_Data_Struct_Pugixml : public Config_Data_Structure
{
public:
	Config_Data_Struct_Pugixml();
	/// Create a new object with a XML tree from a file
	/**
	 * \param file_name name and path of the XML file to parse
	 */
	explicit Config_Data_Struct_Pugixml(::std::string file_name);

	/// Create a new object with a XML tree form a iostream
	/**
	 * \param stream contains the XML data
	 */
	explicit Config_Data_Struct_Pugixml(::std::iostream& stream);

	//
	// Copy constructor (same limitation than the copy operator)
	/*
	 * \param src Object to take for the object creation. The limition from the copy operator is also valid for the copy constructor.
	 */
	//Config_Data_Struct_Pugixml(const Config_Data_Struct_Pugixml& src);


	virtual ~Config_Data_Struct_Pugixml();

	/// Parse a XML file
	void Parse_File(::std::string file_name);

	/// Parse a given stream
	void Parse_Stream(::std::iostream& stream);

	/// Set the current tree pointer to the root of the tree
	void Goto_Root();

	/// Set the current tree pointer to tree with the name param_name
	int Goto_Tree(std::string tree_name, bool error_out = true);

	/// Set the current tree pointer to a subtree with a name
	int Goto_Sub_Tree(std::string tree_name, bool error_out = true);

	/// Read and convert a value with a given name from the current subtree
	/**
	 * \param branch_name name of the branch to read
	 * \param value variable to store the value
	 * \param error_out if true an error is given out if the branch do not exists or exists multiple times
	 * \return 0 if one tag was found, -1 if no tag was found, in case of multiple tags it returns
	 * the number of tags found.
	 * The function throws the invalid_argument exception from stdexcept header file if data in the
	 * tag cannot be converted correctly to the data type. The conversion is very strict,
	 * only if all characters except withespace are used for the conversion, the function does not
	 * throw the exception!
	 */
	int Read_Value(std::string branch_name, int &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, unsigned int &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, long long int &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, unsigned long long int &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, float &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, double &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, std::string &value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, Buffer<int>& value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, Buffer<unsigned int>& value, bool error_out = true);

	/// \overload
	int Read_Value(std::string branch_name, Buffer<float>& value, bool error_out = true);


	/// The value is converted into upper cases.
	/** The read string is converted into upper cases!
	 * \overload
	 */
	int Read_Value(std::string branch_name, bool &value, bool error_out=true);

	/// Read the raw data into a string (useful for e.g., filenames)
	/** same syntax than Read_Value function
	 *
	 */
	int Read_Raw_Value(std::string branch_name, std::string &value, bool error_out=true);

	/// Read the attribute from a given branch from the current subtree.
	/**
	 * \param branch_name name of the branch to read the attribute
	 * \param attrib_name name of the attribute to read
	 * \param value variable to store the attribute value
	 * \param error_out if true an error is given out if the branch or the attribute
	 * do not exists or exists multiple times
	 * \return -1 if the requested attribute does not exist, 0 if no error
	 */
	int Read_Attribute(std::string branch_name, std::string attrib_name, std::string &value, bool error_out = true);


	/// Write the given branch name with the given value into the current subtree.
	/**
	 * All parameter where taken as they are without any conversion!
	 * \param branch_name name of the branch to create
	 * \param value node value of the branch
	 */
	void Write_Raw_Value(const ::std::string& branch_name, const ::std::string& value);

	/// Returns list which contains all name of the subtag of the first tree level.
	Generic_List<std::string, std::vector> Get_All_Subtag_Name();

	/// Streaming operator to print out the current XML tree
	inline friend std::ostream & operator<<(std::ostream & os, Config_Data_Struct_Pugixml& data)
	{
		std::stringstream xml_tree;
		data.doc_.print(xml_tree);
		os << xml_tree.str();
		return os;
	}

	/// Return the modified flag
	bool modified()
	{
		return modified_;
	}

	/// Return the document root object
	const pugi::xml_document& doc()
	{
		return doc_;
	}

	/// Copy operator, copy the internal stored XML tree, but not the current working position in the tree.
	/*Config_Data_Struct_Pugixml &operator=(const Config_Data_Struct_Pugixml& src)
	{
		this->doc_.reset(src.doc_);
		this->modified_ = src.modified_;
		this->current_node_ = this->doc_;

		return *this;
	}*/

protected:

	/// Variable to hold the xml document
	pugi::xml_document doc_;

	/// Current node for processing
	pugi::xml_node current_node_;

	/// Flag to indicate if the tree was changed with a write function
	bool modified_;

	/// Check parse result an make error handling
	void check_parse_result(pugi::xml_parse_result);

	/// Parse the text node of one tag
	int parse_tag_content(pugi::xml_node node, ::std::string tag_name, ::std::string& string_out, bool status_out = true) const;

	/// Check for multiple node with the same name
	int check_multiple_tags(pugi::xpath_node_set node_set, ::std::string tag_name, bool status_out = true) const;

	/// Search for a given node in the subtree of the given node and returns a list of nodes
	pugi::xpath_node_set search_node(pugi::xml_node node, ::std::string tag_name, bool status_out = true) const;

	/// Generates a uniform exception (in case of conversion error)
	void throw_conversion_error(std::exception& e, const std::string& branch_name);
};
#endif /* CONFIG_DATA_STRUCT_PUGIXML_H_ */

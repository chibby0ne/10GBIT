//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the abstract data structure class to hold the configuration tree of one module.
/// \author Timo Lehnigk-Emden
/// \date   2010/11/05
//

#ifndef CONFIG_DATA_STRUCT_H_
#define CONFIG_DATA_STRUCT_H_

#include <string>
#include <stdexcept>
#include "buffer.h"

/// Abstract interface class to hold the configuration tree for one module.
/**
 * This is only an abstract class and cannot be instantiated for that reason.
 * Go to the derived classes to get detailed information for the usage.
 */
class Config_Data_Structure
{
public:
	//Config_Data_Structure() {}
	virtual ~Config_Data_Structure() {}

	/// Set the current tree pointer to the root of the tree
	virtual void Goto_Root() = 0;

	/// Set the current tree pointer to tree with the name param_name
	virtual int Goto_Tree(::std::string tree_name, bool error_out = true) = 0;

	/// Set the current tree pointer to a subtree with a name
	virtual int Goto_Sub_Tree(::std::string tree_name, bool error_out = true) = 0;

	/// Read a value with a given name from the current subtree
	/**
	 * \param branch_name name of the branch to read
	 * \param value variable to store the value
	 * \param error_out if true an error is given out if the branch do not exists or exists multiple times
	 */
	virtual int Read_Value(::std::string branch_name,int &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name,unsigned int &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name,long long int &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name,unsigned long long int &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name,float &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name,double &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name,::std::string &value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name, Buffer<int>& value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name, Buffer<unsigned int>& value, bool error_out = true) = 0;

	/// \overload
	virtual int Read_Value(::std::string branch_name, Buffer<float>& value, bool error_out = true) = 0;

	/// The value is converted into upper cases.
	/**
	 * The read string is converted into upper case!
	 * \overload
	 */
	virtual int Read_Value(::std::string branch_name,bool &value, bool error_out = true) = 0;

	/// Read the raw data into a string (useful for e.g., filenames)
	/**
	 * Same syntax as Read_Value function, return the raw string without conversion.
	 */
	virtual int Read_Raw_Value(::std::string branch_name,::std::string &value, bool error_out=true)=0;

	/// Read the attribute from a given branch from the current subtree
	/**
	 * \param branch_name name of the branch to read the attribute
	 * \param attrib_name name of the attribute to read
	 * \param value variable to store the attribute value
	 * \param error_out if true an error is given out if the branch or the attribute
	 * do not exists or exists multiple times
	 * \return -1 if the requested attribute does not exist, 0 if no error
	 */
	virtual int Read_Attribute(::std::string branch_name, ::std::string attrib_name,::std::string &value, bool error_out=true)=0;

	/// Template function with no functionality to define a needed interface
	template<class T> int Read_Value(::std::string branch_name, T &value, bool error_out=true)
	{
		throw std::runtime_error("Programming error. Function \"Read_Value\" as template may no be used!");
		return -1;
	}

	/// Write the given branch name with the given value into the data structure object at the current position.
	virtual void Write_Raw_Value(const ::std::string& branch_name, const ::std::string& value) = 0;
};
#endif

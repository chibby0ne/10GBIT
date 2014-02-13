//
//  Copyright © 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems
//
/// \file
/// \brief  Definition of the abstract data structure class to hold the result tree
/// \author Timo Lehnigk-Emden
/// \date   2011/01/13
//

#ifndef RESULT_DATA_STRUCT_H_
#define RESULT_DATA_STRUCT_H_

#include <string>

/// Abstract interface class for the result tree
/**
 * \ingroup infrastructure
 */
class Result_Data_Struct
{
public:
	Result_Data_Struct() {}
	virtual ~Result_Data_Struct() {}

	/// Set the current tree pointer to the root of the tree
	virtual void Goto_Root() = 0;

	/// Set the current tree pointer to tree with the name param_name
	virtual int Goto_Tree(::std::string tree_name, bool error_out=true) = 0;

	/// Set the current tree pointer to a subtree with a name
	virtual int Goto_Sub_Tree(::std::string tree_name, bool error_out=true) = 0;

	/// Go up in the tree
	/**
	 * \param levels number of level to go up (default = 1)
	 */
	virtual void Goto_Up(unsigned int levels = 1) = 0;

	/// Write new data
	virtual void Create_Value(::std::string branch_name, ::std::string value) = 0;

	/// Write attribute to the current value
	virtual void Create_Attribute(::std::string attrib_name, ::std::string value) = 0;
};

#endif

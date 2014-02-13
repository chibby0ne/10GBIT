//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of a list to store the pointers of the parameter objects of a module
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef PARAM_LIST_H_
#define PARAM_LIST_H_

#include "../assistance/generic_pointer_list.h"
#include "param_abstract.h"
#include <vector>

namespace cse_lib {

/// Defines a list to store instances of type Param.
/**
 * This list is used print out all Parameter of a module.  In addition, a
 * member variable indicates if a stored instance value has change.
 *
 * \ingroup base_class
 */
class Param_List : public Generic_Pointer_List<Param_Abstract, vector>
{

public:

	Param_List()
	{
		config_modified_ = true;
	}

	virtual ~Param_List(){};

	/// Search for parameter by its name, the first object with the given name is returned.
	/**
	 * \param name of the parameter object pointer to find.
	 * \return pointer to the object with the given name. If no object with the name is in the list
	 *  the zero pointer is returned.
	 */
	Param_Abstract* Search(const string& name)
	{
		Param_Abstract* ptr = 0;

		for(iter_ = obj_list_.begin(); iter_ != obj_list_.end(); iter_++)
		{
			if((*iter_)->name() == name)
			{
				ptr = (*iter_);
				break;
			}
		}

		return ptr;
	}

	/// Set the config_modified_ member.
	void config_modified(bool value)
	{
		config_modified_ = value;
	}

	/// Return a reference to the config_modified_ member.
	/**
	 * The member can be modified outside the class to implement a method to detect the changing
	 * of a parameter of type Param with a single variable.
	 */
	bool& config_modified()
	{
		return config_modified_;
	}

protected:

	/// This flag indicates if the configuration is modified.
	bool config_modified_;
};
}
#endif

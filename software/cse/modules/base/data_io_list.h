//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Defintion of a list to store data I/O object pointers
/// \author Timo Lehnigk-Emden
/// \date   2010/11/05
//
#ifndef DATA_IO_LIST_H_
#define DATA_IO_LIST_H_

#include <vector>
#include "../assistance/generic_pointer_list.h"
#include "data_io_abstract.h"

namespace cse_lib {

/// List to store the input or output data ports of a module
/**
 * Each module has a separate list of the input and output data ports.
 * These lists are necessary for the automatic connection of the ports.
 *
 * \ingroup base_class
 */
class Data_IO_List : public Generic_Pointer_List<Data_IO_Abstract, vector>
{

public:

	Data_IO_List()
	{
		port_modified_ = true;
	}
	virtual ~Data_IO_List(){}


	/// Return the pointer of the buffer object with the given number
	/**
	 * \param id Number of the object in the list to return
	 * \return A pointer to the object is returned, if the requested number is greater than
	 *         the list size the zero pointer is returned.
	 */
	Data_IO_Abstract* Get_Element_By_Id(unsigned int id)
	{
		// Return the zero pointer if no object with this number exists
		if(id >= length())
			return 0;
		else
			return (*(obj_list_.begin() + id));
	}


	/// Return the pointer of the buffer object with the given name of the port
	/**
	 * \param name Name of the object in the list to return
	 * \return A pointer to the object is returned, if the requested object with the name
	 *         does not exits the zero pointer is returned.
	 */
	Data_IO_Abstract* Get_Element_By_Name(string name)
	{
		// Search for the port element with the instance name
		for(iter_ = obj_list_.begin(); iter_ != obj_list_.end(); iter_++)
		{
			if((*iter_)->instance_name() == name)
				return (*iter_);
		}

		// Return the zero pointer if no object with the given name was found
		return 0;
	}

	/// Set the port_modified_ member.
	void port_modified(bool value)
	{
		port_modified_ = value;
	}

	/// Get the port_modified_ member.
	bool& port_modified()
	{
		return port_modified_;
	}

private:

	bool port_modified_;

};
}
#endif

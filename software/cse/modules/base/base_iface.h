//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Abstract base Interface class for data I/O of the modules
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef BASE_INTERF_H_
#define BASE_INTERF_H_

#include "cse_namespace.h"
#include "data_io_list.h"
#include "data_io.h"
#include "status_out_list.h"
#include "status_out.h"
#include "../../assistance/unique_id_instance.h"

namespace cse_lib {

/// Abstract Base Interface Class
/**
 * This class is the parent of all interface classes and therefore of all
 * modules of the simulation chain. The public Run function is mandatory for
 * the scheduler. The Buffer_Event class can be overwritten for i.e. automatic
 * block length support.
 *
 * \ingroup base_class
 */
class Base_Interface
{
public:
	virtual ~Base_Interface() {}

	/// Interface for scheduler, this function processes one block
	virtual int Run() = 0;

	/// Optional function, which can be overloaded.
	/** Modules can overload this functions to keep his values for the next calling of the
	 *  run function. This can be used in iterative systems.
	 *  If the function is called before the Run function, the module will keep its internal state
	 *  for the next call of the Run function.
	 */
	virtual void Keep_Previous_State(){}

	/// Return a reference to the list which contains all status_out objects
	Status_Out_List& status_out_list()
	{
		return status_out_list_;
	}

protected:

	/// List to hold all input data ports
	Data_IO_List input_data_list_;

	/// List to hold all output data ports
	Data_IO_List output_data_list_;

	/// List to hold all status out ports
	Status_Out_List status_out_list_;
};
}
#endif

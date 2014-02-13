//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the abstract data I/O classes of the modules
/// \author Timo Lehnigk-Emden
/// \date   2010/11/05
//
#ifndef DATA_IO_ABSTRACT_H_
#define DATA_IO_ABSTRACT_H_

#include "../assistance/buffer_abstract.h"

namespace cse_lib {
/// Abstract class for the data port input and output of the modules
/**
 * \ingroup base_class
 */
class Data_IO_Abstract
{

public:

	virtual ~Data_IO_Abstract() {}

	/// Abstract function which is overload in the derived classes
	virtual void Set(Buffer_Abstract* in_buffer)=0;

	/// Abstract function which is overload in the derived classes
	virtual Buffer_Abstract* Get() = 0;

	/// Set the name of the port instance
	void instance_name(string instance_name) {instance_name_ = instance_name;}

	/// Return the name of the port instance
	string instance_name() {return instance_name_;}

protected:
	/// Name of the data IO port
	string instance_name_;
};
}
#endif

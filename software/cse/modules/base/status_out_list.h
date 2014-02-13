//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  List definition for the status_out class
/// \author Timo Lehnigk-Emden
/// \date   2011/01/13
//

#ifndef STATUS_OUT_LIST_H_
#define STATUS_OUT_LIST_H_

#include "../assistance/generic_pointer_list.h"
#include "status_out_abstract.h"
#include <vector>

namespace cse_lib {

/// Defines a list to store instances of type status out.
/**
 * This list is used print out all Parameter of a module.  In addition, a
 * member variable indicates if a stored instance value has change.
 *
 * \ingroup base_class
 */
class Status_Out_List : public Generic_Pointer_List<Status_Out_Abstract, vector>
{
public:
	Status_Out_List() {}

	virtual ~Status_Out_List() {}
};
}

#endif

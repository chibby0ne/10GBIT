//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Interface class of the table-based interleaver
/// \author Timo Lehnigk-Emden
/// \date   2011/02/24
//

#ifndef INTERLEAVER_IFACE_H_
#define INTERLEAVER_IFACE_H_


#include "../base/base_iface.h"


namespace cse_lib {

/// Data interface class of the general (de-)interleaver
/**
 * \ingroup interface
 */
template<class T>
class Interleaver_Interface : public Base_Interface
{
public:

	Interleaver_Interface()
	{
		input.Register("input", input_data_list_);
		output.Register("output", output_data_list_);
	}

	virtual ~Interleaver_Interface() { }

	/// Input holding values to (de-)interleave.
	Data_In<T> input;

	/// Output giving (de-)interleaved values.
	Data_Out<T> output;
};
}
#endif

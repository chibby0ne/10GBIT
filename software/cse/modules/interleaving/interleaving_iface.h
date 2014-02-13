//
//  Copyright (C) 2011 creonic UG (haftungsbeschränkt)
//
//  This file is part of the creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Interface of the table based interleaver
/// \author Timo Lehnigk-Emden
/// \date   2011/02/24
//

#ifndef INTERLEAVING_IFACE_H_
#define INTERLEAVING_IFACE_H_


#include "../base/base_iface.h"


namespace cse_lib {

/// Data interface class of the QPP interleaver
/**
 * \ingroup interface
 */
template<class T>
class Interleaving_Interface : public Base_Interface
{
public:

	Interleaving_Interface()
	{
		input.Register("input", input_data_list_);
		output.Register("output", output_data_list_);
	}

	virtual ~Interleaving_Interface() { }

	/// Input holding values to (de-)interleave.
	Data_In<T> input;

	/// Output giving (de-)interleaved values.
	Data_Out<T> output;
};
}
#endif

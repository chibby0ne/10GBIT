//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Encoder interface
/// \author Matthias Alles
/// \date   2010/12/02
//
#ifndef ENCODER_IFACE_H_
#define ENCODER_IFACE_H_

#include "../base/base_iface.h"


namespace cse_lib {

/// Interface class for channel encoders
/**
 * \ingroup interfaces
 */
class Encoder_Interface : public Base_Interface
{

public:

	Encoder_Interface()
	{
		input_bits.Register("input_bits", input_data_list_);
		output_bits.Register("output_bits", output_data_list_);
	}

	virtual ~Encoder_Interface() {}

	Data_In<unsigned int> input_bits;
	Data_Out<unsigned int> output_bits;
};
}
#endif // ENCODER_IFACE_H_

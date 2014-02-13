//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the data interface class of the bit source
/// \author Timo Lehnigk-Emden
/// \date   2010/11/05
//
#ifndef SOURCE_BITS_INTERFACE_H_
#define SOURCE_BITS_INTERFACE_H_

#include "../base/base_iface.h"


namespace cse_lib {

/// Data interface class of the source bits class
/**
 * \ingroup interface
 */
class Source_Bits_Interface : public Base_Interface
{

public:

	Source_Bits_Interface()
	{
		output_bits.Register("output_bits", output_data_list_);
	}

	virtual ~Source_Bits_Interface() {}

	/// Output port which holds the generated bit sequence.
	Data_Out<unsigned int> output_bits;
};
}
#endif

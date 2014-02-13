//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the interface for the float to fixpoint converter module.
/// \author Matthias Alles
/// \date   2010/12/13
//

#ifndef CONV_FLOAT_FIXP_IFACE_H_
#define CONV_FLOAT_FIXP_IFACE_H_

#include "../base/base_iface.h"


namespace cse_lib {

/// Data interface class of the float to fixpoint converter.
/**
 * \ingroup interface
 */
template<class T_in, class T_out> class Converter_Float_Fixpoint_Interface : public Base_Interface
{

public:

	Converter_Float_Fixpoint_Interface()
	{
		input.Register("input", input_data_list_);
		output.Register("output", output_data_list_);
	}

	virtual ~Converter_Float_Fixpoint_Interface() { }

	/// Input holding floating point values.
	Data_In<T_in> input;

	/// Output giving fix-point number representation.
	Data_Out<T_out> output;
};
}
#endif // CONV_FLOAT_FIXP_IFACE_H_

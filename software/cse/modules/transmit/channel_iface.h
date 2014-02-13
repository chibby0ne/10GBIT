//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  SISO channel interface
/// \author Matthias Alles
/// \date   2010/12/01
//

#ifndef CHANNEL_INTERF_H_
#define CHANNEL_INTERF_H_

#include <complex>
#include "../base/base_iface.h"

namespace cse_lib {

/// SISO channel interface class
/**
 * \ingroup interface
 */
class Channel_Interface : public Base_Interface
{

public:

	Channel_Interface()
	{
		input_symb.Register("input_symb", input_data_list_);
		output_symb.Register("output_symb", output_data_list_);
	}

	virtual ~Channel_Interface() {}

	/// Complex input symbols
	Data_In<complex<float> > input_symb;

	/// Complex output symbols with added noise
	Data_Out<complex<float> > output_symb;

};
}
#endif // CHANNEL_INTERF_H_

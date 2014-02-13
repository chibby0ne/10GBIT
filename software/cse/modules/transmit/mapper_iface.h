//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  SISO mapper interface.
/// \author Matthias Alles
/// \date   2010/12/10
//

#ifndef MAPPER_IFACE_H_
#define MAPPER_IFACE_H_

#include "../base/base_iface.h"


namespace cse_lib {

/// Abstract single input single output (SISO) interface class.
/**
 * This abstract class defines the interface for a SISO mapper. The input
 * expects bits. The result are the mapping symbols.
 *
 * \ingroup interface
 */
class Mapper_Interface: public Base_Interface
{

public:

	Mapper_Interface()
	{
		input_bits.Register("input_bits", input_data_list_);
		output_symb.Register("output_symb", output_data_list_);
	}

	virtual ~Mapper_Interface(){}


	/// Input buffer are bits only.
	Data_In<unsigned int> input_bits;

	/// Symbols are given as output.
	Data_Out<complex<float> > output_symb;

};
}
#endif // MAPPER_IFACE_H_

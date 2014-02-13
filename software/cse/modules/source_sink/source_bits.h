//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the binary source module
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef SOURCE_BITS_H_
#define SOURCE_BITS_H_

#include <iostream>

#include "source_bits_iface.h"
#include "source_bits_param.h"


namespace cse_lib {

/// Binary source class
/**
 * This class generates random bit sequences.<br>
 * In addition, the two modes ALL_ZERO and ALL_ONE are supported to get an all
 * zero or all one sequence to perform, e.g., zero code word simulations. In
 * addition, a CRC code can be included in the sequence if the parameter
 * crc_poly is set.
 * \ingroup modules
 */
class Source_Bits: public Source_Bits_Interface,
                   public Source_Bits_Parameter
{

public:

	/// Creates an new instance
	Source_Bits();

	/// Destructor
	virtual ~Source_Bits();

	/// Calculates one block of bits and store them into the output buffer
	int Run();


protected:

	/// initialization of algorithm parameters
	void Init();

	/// Memory allocation
	void Resize_Mem();


private:

	// Generate the random bit sequence
	void Generate_Bits();

	// Calculate CRC bits
	void Add_CRC_Bits();

	// Intern used variables only
	unsigned int crc_length_;
	unsigned int curr_seed_;
};
}
#endif

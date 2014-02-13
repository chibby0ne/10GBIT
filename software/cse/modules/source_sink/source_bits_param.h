//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the parameter of the source bits module
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//

#ifndef SOURCE_BITS_PARAM_H_
#define SOURCE_BITS_PARAM_H_

#include "../base/base_param.h"


namespace cse_lib {

/// Parameter Class for the binary source module
/**
 * In this class all parameters for the binary source are defined. This class
 * should not be instantiated itself, it should be used for derivation of a
 * module.
 * \ingroup parameter
 */
class Source_Bits_Parameter : public Base_Parameter
{

public:

	Source_Bits_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}
	virtual ~Source_Bits_Parameter() {}

	/// Unique identification string of the module.
	static string Unique_ID(){ return "Source_Bits";}


	/**************
	 * Parameters *
	 **************/

	/// Start seed of the random generator for the bit sequence.
	Param<unsigned int> start_seed;

	/// Number of bits to generate.
	Param<unsigned int> num_bits;

	/// CRC code to add to the sequence. Only active if unequal to zero.
	/**
	 * The last bits of the codeword are replaced by the CRC code of the previous bit sequence.
	 * The given polynomial determines the number of bits to replace.
	 * The length of the output sequence is not affected by the CRC feature.
	 * The LSB corresponds to the input flip-flop of the shift register.
	 */
	Param<unsigned int> crc_poly;

	/// Enum to choose the operation mode of the binary source.
	enum MODE_ENUM {RANDOM, ALL_ZERO, ALL_ONE};
	/*! \var RANDOM
	  * The source generates a random sequence of zeros and ones (bit sequence).
	  * \var ALL_ZERO
	  * The source generates a sequence of zeros (e.g., for zero code word simulations).
	  * \var ALL_ONE
	  * The source generates a sequence of ones.
	*/

	/// Operation mode (random, all zero or all one sequence)
	Param<MODE_ENUM> mode;

protected:


	/******************
	 * Default values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - start_seed  : 7878089
	 *  - num_bits    : 1000
	 *  - crc_poly    : 0
	 *  - mode        : RANDOM
	 */
	void Set_Default_Values()
	{
		start_seed.Init(7878089, "start_seed", param_list_);
		num_bits.Init(1000, "num_bits", param_list_);
		crc_poly.Init(0x0, "crc_poly", param_list_);
		mode.Init(RANDOM, "mode", param_list_);
		mode.Link_Value_String(RANDOM, "RANDOM");
		mode.Link_Value_String(ALL_ZERO, "ALL_ZERO");
		mode.Link_Value_String(ALL_ONE, "ALL_ONE");
	}
};
}
#endif

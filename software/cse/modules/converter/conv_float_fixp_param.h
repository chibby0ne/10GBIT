//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the parameter of the float to fixpoint parameters.
/// \author Matthias Alles
/// \date   2010/12/13
//

#ifndef CONV_FLOAT_FIXP_PARAM_H_
#define CONV_FLOAT_FIXP_PARAM_H_

#include "../assistance/config_data_struct.h"
#include "../base/base_param.h"
#include <iostream>


namespace cse_lib {

/// Parameter Class for the float to fixpoint converter module.
/**
 * \ingroup parameter
 */
class Converter_Float_Fixpoint_Parameter : public Base_Parameter
{

public:

	Converter_Float_Fixpoint_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}
	virtual ~Converter_Float_Fixpoint_Parameter() { }


	/// Unique identification string of the module.
	static string Unique_ID(){ return "Converter_Float_Fixpoint";}



	/**************
	 * Parameters *
	 **************/

	/// Number of fractional bits of output value (e.g., 2).
	Param<unsigned int> bw_output_fract;

	/// Bits for representing the output value (e.g., 6).
	Param<unsigned int> bw_output;


protected:

	/******************
	 * Default Values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - bw_output       : 6
	 *  - bw_output_fract : 2
	 */
	void Set_Default_Values()
	{
		bw_output.Init(6, "bw_output", param_list_);
		bw_output_fract.Init(2, "bw_output_fract", param_list_);
	}
};
}
#endif // CONV_FLOAT_FIXP_PARAM_H_

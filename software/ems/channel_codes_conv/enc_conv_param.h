//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * enc_conv_param.h
 *
 *  Created on: 15.02.2011
 *      Author: gimmler
 */

#ifndef ENC_CONV_PARAM_H_
#define ENC_CONV_PARAM_H_

#include "../modules/base/base_param.h"

using namespace cse_lib;

///parameter class for NSC convolutional encoder
class Encoder_Convolutional_Parameter: public Base_Parameter {
public:
	Encoder_Convolutional_Parameter() {
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Encoder_Convolutional_Parameter() {
	}

	///unique name: Encoder_Convolutional
	static string Unique_ID() {
		return "Encoder_Convolutional";
	}

	///forward polynomial for parity1
	Param<unsigned int> forward_polynomial1;
	///forward polynomial for parity2
	Param<unsigned int> forward_polynomial2;

protected:
	///set all parameters to default values
	/**
	 * Default values:
	 *  - forward_polynomial1			: 133
	 *  - forward_polynomial2			: 171
	 */
	virtual void Set_Default_Values() {
		forward_polynomial1.Init(133, "forward_polynomial1", param_list_);
		forward_polynomial2.Init(171, "forward_polynomial2", param_list_);
	}

};

#endif /* ENC_CONV_PARAM_H_ */

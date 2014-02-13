//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * dec_conv_param.h
 *
 *  Created on: 16.02.2011
 *      Author: gimmler
 */

#ifndef DEC_CONV_PARAM_H_
#define DEC_CONV_PARAM_H_

#include "../modules/base/base_param.h"

using namespace cse_lib;

///parameter class for NSC convolutional decoder
class Decoder_Convolutional_Parameter: public Base_Parameter {
public:
	Decoder_Convolutional_Parameter() {
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Decoder_Convolutional_Parameter() {
	}

	///unique name: Decoder_Convolutional
	static string Unique_ID() {
		return "Decoder_Convolutional";
	}

	///forward polynomial for parity1
	Param<unsigned int> forward_polynomial1;
	///forward polynomial for parity2
	Param<unsigned int> forward_polynomial2;
	///check after decoding if output is a valid codeword?
	Param<bool> valid_codeword_check;

	///clipping value for extrinsic output LLRs
	Param<double> output_clipping_value;

	///choice of decoding algorithms
	enum DECODING_ALGORITHM {
		BCJR_MAX_LOG_MAP = 0, BCJR_LOG_MAP
	};
	///decoding algorithm
	Param<DECODING_ALGORITHM> dec_algo;

protected:
	///set all parameters to default values
	/**
	 * Default values:
	 *  - forward_polynomial1			: 133
	 *  - forward_polynomial2			: 171
	 *  - valid_codeword_check			: false
	 *  - output_clipping_value			: 2.0
	 *  - dec_algo 						: BCJR_MAX_LOG_MAP
	 */
	virtual void Set_Default_Values() {
		forward_polynomial1.Init(133, "forward_polynomial1", param_list_);
		forward_polynomial2.Init(171, "forward_polynomial2", param_list_);
		valid_codeword_check.Init(false, "valid_codeword_check", param_list_);
		output_clipping_value.Init(2.0, "output_clipping_value", param_list_);

		dec_algo.Init(BCJR_MAX_LOG_MAP, "dec_algo", param_list_);
		dec_algo.Link_Value_String(BCJR_MAX_LOG_MAP, "BCJR_MAX_LOG_MAP");
		dec_algo.Link_Value_String(BCJR_LOG_MAP, "BCJR_LOG_MAP");
	}

};

#endif /* DEC_CONV_PARAM_H_ */

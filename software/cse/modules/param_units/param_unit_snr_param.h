//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of SNR calculations parameter
/// \author Timo Lehnigk-Emden
/// \date   2010/12/15
//

#ifndef PARAM_UNIT_SNR_PARAM_H_
#define PARAM_UNIT_SNR_PARAM_H_

#include "../base/cse_namespace.h"
#include "../base/base_param.h"

namespace cse_lib{

/// Parameter class for the default SNR / noise variance calculation unit
/**
 * \ingroup parameter
 */
class Param_Unit_SNR_Parameter : public Base_Parameter
{

public:
	Param_Unit_SNR_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Param_Unit_SNR_Parameter() {}

	/// Unique identification string of the module.
	static string Unique_ID(){return "Param_Unit_SNR";}


	/********************
	 * Input Parameters *
	 ********************/

	/// Number of information symbols or bits to take for the variance calculation
	Param<unsigned int> num_info_bits;

	/// Number of all other symbols or bits to transmit to calculate the variance from a given Eb over N0
	Param<unsigned int> num_other_bits;

	/// Number of bits per modulation symbol
	/**
	 * If set to zero, BPSK is assumed with all noise energy on the real part.
	 */
	Param<unsigned int> num_bits_per_symbol;

	/// Input signal to noise ratio (meaning depends on the input mode type)
	Param<float> input_snr;

	/// Enum type to set the input type (Es or Eb over N0)
	enum INPUT_TYPE_ENUM {ES_N0, EB_N0};
	/** \var ES_N0
	 *  Input is transmit signal to noise ratio Es (including all bits to transmit)
	 *  \var EB_N0
	 *  Input is information bits signal to noise ratio Eb
	 */

	/// Interpretation of the input signal to noise ratio
	Param<INPUT_TYPE_ENUM> input_type;


	/*********************
	 * Output Parameters *
	 *********************/

	/// Output noise variance
	Param_Out<float> noise_variance;

	/// Output Eb over N0 to calculate the noise variance for a complex channel
	Param_Out<float> eb_n0;

	/// Output Signal to noise ratio in dB (Es over N0)
	Param_Out<float> es_n0;


protected:


	/******************
	 * Default values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - input_snr           : 0.0
	 *  - input_type          : EB_N0
	 *  - num_info_bits       : 1
	 *  - num_other_bits      : 1
	 *  - num_bits_per_symbol : 1
	 */
	void Set_Default_Values()
	{
		input_snr.Init(0.0, "input_snr", param_list_);
		input_type.Init(EB_N0, "input_type", param_list_);
		num_info_bits.Init(1, "num_info_bits", param_list_);
		num_other_bits.Init(1, "num_other_bits", param_list_);
		num_bits_per_symbol.Init(1, "num_bits_per_symbol", param_list_);

		// Register output parameter in the parameter list
		// String is to link the parameter with the XML tags!
		noise_variance.Register("noise_variance", param_out_list_);
		eb_n0.Register("eb_n0", param_out_list_);
		es_n0.Register("es_n0", param_out_list_);

		// Link enum types to string
		input_type.Link_Value_String(ES_N0, "ES_N0");
		input_type.Link_Value_String(EB_N0, "EB_N0");
	}
};

}
#endif

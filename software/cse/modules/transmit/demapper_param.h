//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the parameters for the SISO demapper.
/// \author Timo Lehnigk-Emden
/// \date   2010/12/10
//
#ifndef DEMAPPER_PARAM_H_
#define DEMAPPER_PARAM_H_

#include "../base/base_param.h"
#include "mapper_share.h"

namespace cse_lib {

/// Parameter class for the SISO demapper.
/**
 * \ingroup parameter
 */
class Demapper_Parameter : public Base_Parameter
{

public:

	Demapper_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Demapper_Parameter(){}

	/// Unique identification string of the module.
	static string Unique_ID(){return "Demapper";}


	/**************
	 * Parameters *
	 **************/

	/// Modulation scheme
	Param<Mapper_Share::MAPPINGS_ENUM> mapping;

	/// Enum type for the channel reliability factor behavior
	enum CHANNEL_RELIABILITY_ENUM {OPTIMAL_CRF, FIXED_CRF};
	/*! \var OPTIMAL_CRF
	 *  The optimal noise variance is taken for demapping
	 *  \var FIXED_CRF
	 *  The given channel reliability factor is taken for demappping
	 */

	/// Channel reliability behavior, LLR calculation behavior
	Param<CHANNEL_RELIABILITY_ENUM> channel_reliability;

	/// Noise variance used for demapping, this is the overall variance (imaginary and real part)
	/**
	 * The given noise variance is equally split to the two dimension.
	 * For the real and the imaginary part respectively a noise variance of noise_variance / 2
	 * is taken for the demapping process.
	 */
	Param<float> noise_variance;

	/// Channel reliability factor which is used in case of channel_reliability = FIXED_CRF
	Param<float> channel_reliability_factor;


	/// PSK radius ratio, only relevant for 16-APSK and 32-APSK modulation.
	/**
	 * This ratio gives the radius of the inner circle to the next circle. By
	 * default, 16-APSK uses two circles of constellation points.
	 * E.g., in DVB-S2, it is mandatory to support different ratios for different
	 * coding rates.
	 */
	Param<float> psk_radius_ratio_r2;


	/// PSK radius ratio, only relevant for 32-APSK modulation.
	/**
	 * This ratio gives the radius of the inner circle to the outer circle.  By
	 * default 32-APSK uses three circles of constellation points.  E.g., in
	 * DVB-S2, it is mandatory to support different ratios for different coding
	 * rates.
	 */
	Param<float> psk_radius_ratio_r3;


	/*****************
	 * Configuration *
	 *****************/

	/// Configure the module from the given data structure.
	void Config_From_Data_Struct(Config_Data_Structure& config_tree)
	{
		mapping.Set_By_Data_Struct(config_tree);

		channel_reliability.Set_By_Data_Struct(config_tree);

		// Dependent parameter
		if(channel_reliability() == OPTIMAL_CRF)
			noise_variance.Set_By_Data_Struct(config_tree);
		else if(channel_reliability() == FIXED_CRF)
			channel_reliability_factor.Set_By_Data_Struct(config_tree);

		// Optional parameters

		// 16-APSK or 32-APSK => get radius ratio r2
		if(mapping() == Mapper_Share::MAP_16_APSK || mapping() == Mapper_Share::MAP_32_APSK)
			psk_radius_ratio_r2.Set_By_Data_Struct(config_tree);

		// 32-APSK => get radius ratio r3
		if(mapping() == Mapper_Share::MAP_32_APSK)
			psk_radius_ratio_r3.Set_By_Data_Struct(config_tree);
	}


protected:


	/******************
	 * Default values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - mapping                    : MAP_BPSK
	 *  - channel_reliability        : OPTIMAL_CRF
	 *  - noise_variance             : 1.0
	 *  - channel_reliability_factor : 1.0
	 *  - psk_radius_ratio_r2        : 2.5
	 *  - psk_radius_ratio_r3        : 4.5
	 */
	void Set_Default_Values()
	{
		mapping.Init(Mapper_Share::MAP_BPSK, "mapping", param_list_);
		channel_reliability.Init(OPTIMAL_CRF, "channel_reliability", param_list_);
		noise_variance.Init(1.0, "noise_variance", param_list_);
		channel_reliability_factor.Init(1.0, "channel_reliability_factor", param_list_);
		psk_radius_ratio_r2.Init(2.5, "psk_radius_ratio_r2", param_list_);
		psk_radius_ratio_r3.Init(4.5, "psk_radius_ratio_r3", param_list_);

		// Link enum types to strings
		mapping.Link_Value_String(Mapper_Share::Mapping_Enum_Str());
		channel_reliability.Link_Value_String(OPTIMAL_CRF, "OPTIMAL_CRF");
		channel_reliability.Link_Value_String(FIXED_CRF, "FIXED_CRF");
	}
};
}
#endif

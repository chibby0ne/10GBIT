//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  SISO mapper parameters.
/// \author Matthias Alles
/// \date   2010/12/10
//

#ifndef MAPPER_PARAM_H_
#define MAPPER_PARAM_H_

#include <complex>
#include "../base/base_param.h"
#include "mapper_share.h"


namespace cse_lib {

/// Parameter class for the SISO mapper.
/**
 *
 * \ingroup parameter
 */
class Mapper_Parameter : public Base_Parameter
{

public:

	Mapper_Parameter()
	{
		Set_Default_Values();
		instance_name(Unique_ID());
	}

	virtual ~Mapper_Parameter(){}

	/// Unique identification string of the module.
	static string Unique_ID(){return "Mapper";}


	/**************
	 * Parameters *
	 **************/

	/// Modulation scheme
	Param<Mapper_Share::MAPPINGS_ENUM> mapping;


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

	void Config_From_Data_Struct(Config_Data_Structure& config_tree)
	{
		mapping.Set_By_Data_Struct(config_tree);

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
	 *  - psk_radius_ratio_r2        : 2.5
	 *  - psk_radius_ratio_r3        : 4.5
	 */
	void Set_Default_Values()
	{
		mapping.Init(Mapper_Share::MAP_BPSK, "mapping", param_list_);
		mapping.Link_Value_String(Mapper_Share::Mapping_Enum_Str());
		psk_radius_ratio_r2.Init(2.5, "psk_radius_ratio_r2", param_list_);
		psk_radius_ratio_r3.Init(4.5, "psk_radius_ratio_r3", param_list_);
	}

};
}
#endif // MAPPER_PARAM_H_

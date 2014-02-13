//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * mimo_channel_rayleigh_param.h
 *
 *  Created on: 03.02.2011
 *      Author: gimmler
 */

#ifndef MIMO_CHANNEL_RAYLEIGH_PARAM_H_
#define MIMO_CHANNEL_RAYLEIGH_PARAM_H_

#include "../modules/base/base_param.h"

using namespace cse_lib;

///parameter class for Rayleigh fading mimo channel
class Mimo_Channel_Rayleigh_Parameter: public Base_Parameter {
public:
	Mimo_Channel_Rayleigh_Parameter() {
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Mimo_Channel_Rayleigh_Parameter() {
	}

	///unique name: Mimo_Channel_Rayleigh
	static string Unique_ID() {
		return "Mimo_Channel_Rayleigh";
	}

	///number of transmit/receive antennas
	Param<int> num_antennas;
	///number of channel uses for which the channel matrix stays constant
	Param<int> num_constant_channels;
	///signal-to-noise ratio on symbol-layer
	Param<float> es_n0;
	///transmit with noise? (choose false for debugging)
	Param<bool> noise_on;
	///enum for all fading types
	enum FADING_TYPE {
		QUASI_STATIC = -1, BLOCK_FADING, ERGODIC
	};
	///channel fading
	Param<FADING_TYPE> fading;

protected:

	/**
	 * \brief set all parameters to default values
	 *
	 * Default values:
	 *  - num_antennas				: 4
	 *  - num_constant_channels		: 1
	 *  - es_n0						: 1.0
	 *  - noise_on					: true
	 *  - fading					: ERGODIC
	 */
	virtual void Set_Default_Values() {
		num_antennas.Init(4, "num_antennas", param_list_);
		num_constant_channels.Init(1, "num_constant_channels", param_list_);
		es_n0.Init(1.0, "es_n0", param_list_);
		noise_on.Init(true, "channel_noise_on", param_list_);
		fading.Init(ERGODIC, "fading", param_list_);

		fading.Link_Value_String(QUASI_STATIC, "QUASI_STATIC");
		fading.Link_Value_String(ERGODIC, "ERGODIC");
		fading.Link_Value_String(BLOCK_FADING, "BLOCK_FADING");
	}

};

#endif /* MIMO_CHANNEL_RAYLEIGH_PARAM_H_ */


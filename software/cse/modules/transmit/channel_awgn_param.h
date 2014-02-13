//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  AWGN channel parameters
/// \author Matthias Alles
/// \date   2010/12/01
//

#ifndef CHANNEL_AWGN_PARAM_H_
#define CHANNEL_AWGN_PARAM_H_

#include "../base/base_param.h"

namespace cse_lib {

/// Parameter Class for a SISO AWGN channel
/**
 * \ingroup parameter
 */
class Channel_AWGN_Parameter: public Base_Parameter
{

public:

	Channel_AWGN_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Channel_AWGN_Parameter(){}

	/// Unique identification string of the module.
	static string Unique_ID(){return "Channel_AWGN";}


	/**************
	 * Parameters *
	 **************/

	/// Start seed of the noise generator
	Param<unsigned int> start_seed;

	/// Variance of the channel. The variance is separated to the two dimensions (I and Q).
	Param<float> noise_variance;


protected:


	/******************
	 * Default values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - start_seed     : 12424
	 *  - noise_variance : 1.0
	 */
	void Set_Default_Values()
	{
		start_seed.Init(12424, "start_seed", param_list_);
		noise_variance.Init(1.0, "noise_variance", param_list_);
	}
};
}
#endif // CHANNEL_AWGN_PARAM_H_

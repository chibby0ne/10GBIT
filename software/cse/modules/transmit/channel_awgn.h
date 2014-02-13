//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  SISO AWGN channel module definition
/// \author Matthias Alles
/// \date   2010/12/01
//

#ifndef CHANNEL_AWGN_H_
#define CHANNEL_AWGN_H_

#include <iostream>
#include "channel_awgn_param.h"
#include "channel_iface.h"
#include "../assistance/gaussian_random_generator.h"

namespace cse_lib {

/// SISO AWGN Channel
/**
 * \ingroup modules
 */
class Channel_AWGN: public Channel_Interface,
                    public Channel_AWGN_Parameter
{

public:

	Channel_AWGN();
	virtual ~Channel_AWGN();

	int Run();

private:

	void Init();
	void Add_AWGN();
	// External random generator
	Gaussian_Random_Generator rand_gen_;

	// Std variation, calculated from the variance
	float std_variation_;

};
}
#endif // CHANNEL_AWGN_H_

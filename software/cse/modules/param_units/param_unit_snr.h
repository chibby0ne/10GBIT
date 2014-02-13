//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of SNR calculations
/// \author Timo Lehnigk-Emden
/// \date   2010/12/13
//

#ifndef PARAM_UNIT_SNR_H_
#define PARAM_UNIT_SNR_H_

#include "../base/cse_namespace.h"
#include "../base/base_iface.h"
#include "param_unit_snr_param.h"

namespace cse_lib {

/// Standard SNR to noise variance calculation
/**
 * \ingroup param_unit
 */
class Param_Unit_SNR : public Base_Interface,
                       public Param_Unit_SNR_Parameter
{

public:

	Param_Unit_SNR();
	virtual ~Param_Unit_SNR();

	/// Calculate the Eb/N0, Es/N0 and the corresponding variance with the given parameters.
	int Run();

private:
	// Calculate the EB/N0, ES/NO and the noise variance with the default method
	void Calc_SNR_Variance();
};
}
#endif

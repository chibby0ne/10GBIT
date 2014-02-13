//
//  Copyright (C) 2011 - 2012 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Binary LDPC Encoder for IEEE_802_11ad LDPC codes, method definitions
/// \author Markus Fehrenz
/// \date   2012/01/24
//

#include "enc_ldpc_ieee_802_11ad.h"
#include "enc_ldpc_ieee_802_11ad_codes.hpp"

using namespace hlp_fct::logging;
using namespace cse_lib::ieee_802_11ad_codes;
using namespace std;

namespace cse_lib {

void Encoder_LDPC_IEEE_802_11ad::Init()
{
	// Parameterize Encoder
	Set_LDPC_Code_Parameters();

	// Resize output buffer
	try
	{
		// output memory corresponds to codeword size
		output_bits().Resize((k_b_ + m_b_) * z_);
	}
	catch(bad_alloc&)
	{
		Msg(ERROR, instance_name(), "Memory allocation failure!");
		throw;
	}

	param_list_.config_modified(false);
	input_data_list_.port_modified(false);
}


/// Set the parameters for the Run() function of the parent class.
void Encoder_LDPC_IEEE_802_11ad::Set_LDPC_Code_Parameters()
{
	switch(ldpc_code())
	{

	case IEEE_802_11AD_P42_N672_R050:
		k_b_ = 8;
		m_b_ = 8;
		z_   = 42;
		H_   = &H_ieee_802_11ad_p42_n672_r050[0][0];
		break;

	case IEEE_802_11AD_P42_N672_R062:
		k_b_ = 10;
		m_b_ = 6;
		z_   = 42;
		H_   = &H_ieee_802_11ad_p42_n672_r062[0][0];
		break;

	case IEEE_802_11AD_P42_N672_R075:
		k_b_ = 12;
		m_b_ = 4;
		z_   = 42;
		H_   = &H_ieee_802_11ad_p42_n672_r075[0][0];
		break;

	case IEEE_802_11AD_P42_N672_R081:
		k_b_ = 13;
		m_b_ = 3;
		z_   = 42;
		H_   = &H_ieee_802_11ad_p42_n672_r081[0][0];
		break;

	}
};


int Encoder_LDPC_IEEE_802_11ad::Run()
{

	// In case of modified parameter, call local Init() method.
	if(param_list_.config_modified())
		Init();

	Encode_Binary_Lower_Triangular(input_bits(), output_bits());

	return 0;
}
}


//
//  Copyright (C) 2011 - 2012 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Binary LDPC Encoder for IEEE_802_11ad LDPC codes, module definition
/// \author Markus Fehrenz
/// \date   2012/01/24
//

#ifndef ENC_LDPC_IEEE_802_11AD_H_
#define ENC_LDPC_IEEE_802_11AD_H_

#include "enc_iface.h"
#include "enc_ldpc_bin_struct_share.h"
#include "enc_ldpc_ieee_802_11ad_param.h"


namespace cse_lib {

/// LDPC Encoder for IEEE_802_11ad standard (N = 672, P = 24, R = (0.5 | 0.625 | 0.75 | 0.81) )
/**
 * \ingroup modules
 */
class Encoder_LDPC_IEEE_802_11ad : public Encoder_Interface,
                                   public Encoder_LDPC_IEEE_802_11ad_Parameter,
                                   public Encoder_LDPC_Binary_Struct_Share
{

public:

	Encoder_LDPC_IEEE_802_11ad() { };
	virtual ~Encoder_LDPC_IEEE_802_11ad() { };

	int Run();

private:

	void Init();
	void Set_LDPC_Code_Parameters();

};
}
#endif // ENC_LDPC_IEEE_802_11AD_H_

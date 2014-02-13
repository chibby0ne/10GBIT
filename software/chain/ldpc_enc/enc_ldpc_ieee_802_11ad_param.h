//
//  Copyright (C) 2011 - 2012 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Parameter class for IEEE 802.11ad LDPC encoder
/// \author Markus Fehrenz
/// \date   2012/01/24
//

#ifndef ENC_LDPC_IEEE_802_11AD_PARAM_H_
#define ENC_LDPC_IEEE_802_11AD_PARAM_H_

//#include "../base/base_param.h"
#include "cse_lib.h"


namespace cse_lib {

/// Binary LDPC encoder parameter class for IEEE 802.11ad LDPC codes
/**
 * \ingroup parameter
 */
class Encoder_LDPC_IEEE_802_11ad_Parameter : public Base_Parameter
{

public:

	Encoder_LDPC_IEEE_802_11ad_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Encoder_LDPC_IEEE_802_11ad_Parameter(){}

	/// Unique identification string of the module.
	static std::string Unique_ID(){return "Encoder_LDPC_IEEE_802_11ad";}


	/**************
	 * Parameters *
	 **************/

	/// List of available LDPC codes
	enum LDPC_CODE {
		IEEE_802_11AD_P42_N672_R050,
		IEEE_802_11AD_P42_N672_R062,
		IEEE_802_11AD_P42_N672_R075,
		IEEE_802_11AD_P42_N672_R081
	};

	/// Selects LDPC code
	Param<LDPC_CODE> ldpc_code;


protected:


	/******************
	 * Default values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - ldpc_code      : IEEE_802_11AD_P42_N672_R050
	 */
	void Set_Default_Values()
	{
		ldpc_code.Init(IEEE_802_11AD_P42_N672_R050, "ldpc_code", param_list_);
		ldpc_code.Link_Value_String(IEEE_802_11AD_P42_N672_R050, "IEEE_802_11AD_P42_N672_R050");
		ldpc_code.Link_Value_String(IEEE_802_11AD_P42_N672_R062, "IEEE_802_11AD_P42_N672_R062");
		ldpc_code.Link_Value_String(IEEE_802_11AD_P42_N672_R075, "IEEE_802_11AD_P42_N672_R075");
		ldpc_code.Link_Value_String(IEEE_802_11AD_P42_N672_R081, "IEEE_802_11AD_P42_N672_R081");
	}

};
}
#endif // ENC_LDPC_IEEE_802_11_AD_PARAM_H_

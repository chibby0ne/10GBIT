//
//  Copyright (C) 2012 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  IEEE 802.3an LDPC decoder parameters
/// \author Matthias Alles
/// \date   2012/03/26
//

#ifndef DEC_LDPC_IEEE_802_3AN_PARAM_H_
#define DEC_LDPC_IEEE_802_3AN_PARAM_H_

//#include "base/base_param.h"
#include "../../cse/modules/base/base_param.h"
#include "dec_ldpc_bin_hw_share.h"
#include "dec_ldpc_check_node_share.h"


namespace cse_lib {

/// Parameter class for the binary hardware LDPC decoder
/**
 * \ingroup parameter
 */
class Decoder_LDPC_IEEE_802_3an_Parameter : public Base_Parameter
{

public:

	Decoder_LDPC_IEEE_802_3an_Parameter()
	{
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Decoder_LDPC_IEEE_802_3an_Parameter(){}

	/// Unique identification string of the module.
	static std::string Unique_ID(){ return "Decoder_LDPC_IEEE_802_3an";}


	/****************************
	 ** Algorithmic parameters **
	 ****************************/

	/// LDPC decoding scheduling types
	enum DEC_SCHEDULING_ENUM
	{
		TWO_PHASE,         /*!< Two Phase decoding */
		LAYERED            /*!< Layered scheduling (aka Turbo Decoding Message Passing) */
	};

	/// LDPC Decoder scheduling
	Param<DEC_SCHEDULING_ENUM> scheduling;

	/// Number of Decoder Iterations (default: 10)
	Param<unsigned int> num_iterations;

	/// Decoding Algorithm
	Param<Decoder_LDPC_Check_Node_Share::CHECK_NODE_ENUM> dec_algorithm;

	/// Number of messages to take in account for lambda-min algorithm
	Param<unsigned int> num_lambda_min;

	/// Extrinsic scaling factor for Min-Sum algorithm
	Param<float> esf_factor;

	/// Perform parity check on APP values (true) or on extrinsic values (false)
	Param<bool> app_parity_check;

 	/// Number of Partitions for Multi Split-Row
    Param<unsigned int> num_partitions;

    /// Threshold for Split Row
    Param<unsigned int> threshold;

	/******************
	 ** Quantization **
	 ******************/

	/// Number of bits for fractional part of all values (default: 2 bits)
	Param<unsigned int> bw_fract;

	/// Number of bits for all extrinsic values (default: 6 bits)
	Param<unsigned int> bw_extr;

	/// Number of bits for the APP values (default: 9 bits)
	Param<unsigned int> bw_app;


	/*********************
	 ** Code Parameters **
	 *********************/
    enum LDPC_CODE {
        IEEE_802_3AN_P1_N2048_R084
    };

    Param<LDPC_CODE> ldpc_code;

protected:


	/******************
	 * Default Values *
	 ******************/

	/// Set all parameters to their default values.
	/**
	 * Default values:
	 *  - bw_fract         : 2
	 *  - bw_extr          : 6
	 *  - bw_app           : 9
	 *  - num_iterations   : 10
	 *  - num_lambda_min   : 3
	 *  - esf_factor       : 0.875
	 *  - dec_algorithm    : MIN_SUM
	 *  - scheduling       : LAYERED
	 *  - app_parity_check : true
	 */
	void Set_Default_Values()
	{
        ldpc_code.Link_Value_String(IEEE_802_3AN_P1_N2048_R084, "IEEE_802_3AN_P1_N2048_R084");

		dec_algorithm.Link_Value_String(Decoder_LDPC_Check_Node_Share::LAMBDA_MIN, "LAMBDA_MIN");
		dec_algorithm.Link_Value_String(Decoder_LDPC_Check_Node_Share::MIN_SUM,    "MIN_SUM");
		dec_algorithm.Link_Value_String(Decoder_LDPC_Check_Node_Share::MIN_SUM_SELF_CORRECTING, "MIN_SUM_SELF_CORRECTING");
		dec_algorithm.Link_Value_String(Decoder_LDPC_Check_Node_Share::SPLIT_ROW,    "SPLIT_ROW");
		dec_algorithm.Link_Value_String(Decoder_LDPC_Check_Node_Share::SPLIT_ROW_IMPROVED,    "SPLIT_ROW_IMPROVED");
		dec_algorithm.Link_Value_String(Decoder_LDPC_Check_Node_Share::SPLIT_ROW_SELF_CORRECTING,    "SPLIT_ROW_SELF_CORRECTING");
		scheduling.Link_Value_String(TWO_PHASE, "TWO_PHASE");
		scheduling.Link_Value_String(LAYERED,   "LAYERED");

		bw_fract.Init(2, "bw_fract", param_list_);
		bw_extr.Init (6, "bw_extr", param_list_);
		bw_app.Init  (9, "bw_app", param_list_);

		num_iterations.Init(15, "num_iterations", param_list_);
		num_lambda_min.Init(3, "num_lambda_min", param_list_);
		esf_factor.Init(0.875, "esf_factor", param_list_);
		app_parity_check.Init(true, "app_parity_check", param_list_);
 		num_partitions.Init(2, "num_partitions", param_list_);
        threshold.Init(12, "threshold", param_list_);


		dec_algorithm.Init(Decoder_LDPC_Check_Node_Share::MIN_SUM, "dec_algorithm", param_list_);
		scheduling.Init(LAYERED, "scheduling", param_list_);

        ldpc_code.Init(IEEE_802_3AN_P1_N2048_R084, "ldpc_code", param_list_);
	}

};
}
#endif // DEC_LDPC_IEEE_802_3AN_PARAM_H_

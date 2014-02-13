//
//  Copyright (C) 2012 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Check node algorithm definitions
/// \author Matthias Alles
/// \date   2012/11/26
//

#ifndef DEC_LDPC_CHECK_NODE_SHARE_H_
#define DEC_LDPC_CHECK_NODE_SHARE_H_

namespace cse_lib {

/// Parameter class for the binary hardware LDPC decoder
/**
 * \ingroup parameter
 */
class Decoder_LDPC_Check_Node_Share
{

public:

	Decoder_LDPC_Check_Node_Share(){};

	virtual ~Decoder_LDPC_Check_Node_Share(){};

	/// LDPC check node algorithms
	enum CHECK_NODE_ENUM {
		LAMBDA_MIN,      /*!< Lambda-Min, where lambda is given by num_lambda_min */
		MIN_SUM,         /*!< Min-Sum decoding, with extrinsic scaling factor esf_factor */
		LAMBDA_3MIN_3MAG, /*!< Lambda-Min, searching three minima, output three different magnitudes */
        MIN_SUM_SELF_CORRECTING,
        SPLIT_ROW,
        SPLIT_ROW_IMPROVED,
        SPLIT_ROW_SELF_CORRECTING
	};

};
}
#endif

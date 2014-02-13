//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Binary LDPC Encoder for structured codes, shared functionality.
/// \author Matthias Alles
/// \date   2010/12/02
//

#ifndef ENC_LDPC_BIN_STRUCT_H_
#define ENC_LDPC_BIN_STRUCT_H_

#include "enc_iface.h"


namespace cse_lib {

/// Shared functionality for binary LDPC encoder for structured codes.
/**
 * \ingroup share
 */
class Encoder_LDPC_Binary_Struct_Share
{

protected:

	Encoder_LDPC_Binary_Struct_Share() { }
	virtual ~Encoder_LDPC_Binary_Struct_Share() { }


	/***********************************************************
	 * Parameters required by functions of this class to work *
	 ***********************************************************/

	/// Number of submatrices corresponding to systematic information
	u_int k_b_;

	/// Number of submatrices corresponding to parity information
	u_int m_b_;

	/// Submatrix size of parity check matrix
	u_int z_;

	/// Pointer to parity check matrix structure
	int *H_;


	/// If H_ in a standard is defined with left shifts, translate them to right shifts.
	/**
	 * Some standards define the submatrices in the parity check matrix with cyclic
	 * left shifts, some with cyclic right shift. Since right shifts are more common,
	 * the encoder functions of this class use cyclic right shifts. This function
	 * allows for the adaption of cyclic left shifts when defined in the standard to
	 * cyclic right shifts, i.e., it is possible to store the parity check matrix as
	 * defined in the standard, then call this function, which translates the left
	 * shifts to right shifts, and finally call the encoder function.
	 */
	void Change_Shift_Direction(bool lower_triangular_form);


	/// Encoder function for structured LDPC codes.
	/**
	 * Function that performs LDPC encoding of structured LDPC codes as used in
	 * WiMax (802.16e), WiFi (802.11n), WiMedia 1.5, CMMB. The formulas for this
	 * encoding procedure are described in the WiMax standard. Shift values in
	 * H_ are assumed to be cyclic right shifts (adapt H_ when standard
	 * defines right shifts).
	 * The function requires H_ to contain only the information nodes and the
	 * first parity node group.
	 */
	void Encode_Binary_Struct(Buffer<unsigned int> &input_bits,
	                          Buffer<unsigned int> &output_bits);


	/// Encoder function for structured LDPC codes in lower triangular form.
	/**
	 * Function that performs LDPC encoding of structured LDPC codes in lower
	 * triangular form, e.g., used in IEEE 802.15.3c. The encoding of these
	 * codes is straight forward. Shift values in H_ are assumed to be
	 * cyclic right shifts (adapt H_ when standard defines right shifts).
	 * The function requires H_ to be complete.
	 */
	void Encode_Binary_Lower_Triangular(Buffer<unsigned int> &input_bits,
	                                    Buffer<unsigned int> &output_bits);
};
}
#endif // ENC_LDPC_BIN_STRUCT_H_

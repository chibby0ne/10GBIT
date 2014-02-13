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

#include "enc_ldpc_bin_struct_share.h"

namespace cse_lib {


void Encoder_LDPC_Binary_Struct_Share::Change_Shift_Direction(bool lower_triangular_form)
{

	u_int num_vng;        // Number of variable node groups stored in the parity check matrix H_.
	u_int num_cng = m_b_; // Number of check node groups.

	if (lower_triangular_form)
		num_vng = k_b_ + m_b_;
	else
		num_vng = k_b_ + 1;

	for (u_int submatrix = 0; submatrix < num_cng * num_vng; submatrix++)
	{
		int h_entry = *(H_ + submatrix);

		if (h_entry > 0)
			*(H_ + submatrix) = z_ - h_entry;
	}
}


void Encoder_LDPC_Binary_Struct_Share::Encode_Binary_Struct(Buffer<unsigned int> &input_bits,
                                                            Buffer<unsigned int> &output_bits)
{
	Buffer<int, 2> u(k_b_, z_);
	Buffer<int, 2> v(m_b_, z_);
	u_int i, j, k, q;
	int h_entry;

	for (j = 0; j < k_b_; j++)
		for (i = 0; i < z_; i++)
		{
			// init u-vector
			u[j][i] = input_bits[i + z_ * j];
		}

	for (j = 0; j < m_b_; j++)
		for (i = 0; i < z_; i++)
		{
			// init sum_vector
			v[j][i] = 0;
		}

	/*************
	 * Method 1a *
	 *************/

	// calculate v_0
	for (j = 0; j < k_b_; j++)
		for (q = 0; q < m_b_; q++)
		{
			h_entry = *(H_ + q * (k_b_ + 1) + j);

			if( h_entry != -1)
				for(k = 0; k < z_; k++)
					v[0][k] = (v[0][k] + u[j][ (h_entry + k) % z_ ]);
		}

	/*
	 * Warning: it seems that this code does not take into account that the
	 * above equations calculate
	 * P_ p (x, k) v(0)
	 * So, if the positive matrix element in the middle of column k_b is not
	 * zero, but a different positive value, the equations calculate a shifted
	 * version of v(0). One additional reverse shift is necessary for the
	 * correct v(0) calculation.
	 *
	 * WiMax code rate 3/4B will most likely not work with this code.  All
	 * other WiMax codes, WiMedia 1.5 codes, and WiFi (802.11n) codes will
	 * work.
	 */


	// calculate v_1
	for(j = 0; j < k_b_; j++)
	{
		h_entry = *(H_ + 0 * (k_b_ + 1) + j);

		if (h_entry != -1)
		{
			for (k = 0; k < z_; k++) {
				v[1][k] = (v[1][k] + u[j][ (h_entry + k) % z_ ] );
			}
		}
	}

	h_entry = *(H_ + 0 * (k_b_ + 1) + k_b_);

	if (h_entry != -1)
		for (k = 0; k < z_; k++)
			v[1][k] = (v[1][k] + v[0][ (h_entry + k) % z_ ] );

	// calculate v_(i+1)
	for(i = 1; i < m_b_ - 1; i++)
	{
		for(k = 0; k < z_; k++)
		{
			v[i + 1][k] = v[i][k];
		}

		for(j = 0; j < k_b_; j++)
		{
			h_entry = *(H_ + i * (k_b_ + 1) + j);

			if (h_entry != -1)
				for (k = 0; k < z_; k++)
					v[i + 1][k] = (v[i + 1][k] + u[j][ (h_entry + k) % z_ ] );
		}

		h_entry = *(H_ + i * (k_b_ + 1) + k_b_);

		if (h_entry != -1)
			for(k = 0; k < z_; k++)
				v[i + 1][k] = (v[i + 1][k] + v[0][ (h_entry + k) % z_ ] );
	}

	// write the information bits
	for (j = 0; j < k_b_ * z_; j++)
		output_bits[j] = input_bits[j];

	// write the parity bits
	for(i = 0; i < m_b_; i++)
		for(k = 0; k < z_; k++)
		{
			v[i][k] = v[i][k] % 2;
			output_bits[(k_b_ + i) * z_ + k] = v[i][k];
		}

}



void Encoder_LDPC_Binary_Struct_Share::Encode_Binary_Lower_Triangular(Buffer<unsigned int> &input_bits,
	                                                                  Buffer<unsigned int> &output_bits)
{

	u_int num_vngs = k_b_ + m_b_;    // Number of variable node groups.
	int h_entry;

	output_bits.Clear();
	for(u_int i = 0; i < input_bits.length(); i++)
		output_bits[i] = input_bits[i];


	for (u_int cng = 0; cng < m_b_; cng++)
	{
		Buffer<u_int> parity_bits(z_);
		parity_bits.Clear();

		/*
		 * Iterate for each check node group over the corresponding variable
		 * node groups.  For each check node group, one additional variable
		 * node group is required to be taken into account, because of lower
		 * triangular form.
		 */
		for (u_int vng = 0; vng < k_b_ + cng; vng++)
		{
			h_entry = *(H_ + cng * num_vngs + vng);

			// Accumulate parity bits, except for last entry in check node group.
			if (h_entry > -1)
			{
				for (u_int edge_cnt = 0; edge_cnt < z_; edge_cnt++)
				{
					parity_bits[edge_cnt] += output_bits[vng * z_ + (h_entry + edge_cnt) % z_];
				}
			}
		}

		/*
		 * Perform the final assignment of the calculated parity bits to the
		 * output vector.
		 */
		u_int current_vng = k_b_ + cng;    // current variable node group
		u_int cw_bit = current_vng * z_;   // current bit of the codeword we compute

		h_entry = *(H_ + cng * num_vngs + current_vng);

		// Copy back the just calculated parity bits to the output.
		for (u_int edge_cnt = 0; edge_cnt < z_; edge_cnt++)
		{
			output_bits[cw_bit] = parity_bits[(z_ - h_entry + edge_cnt) % z_] % 2;
			cw_bit++;
		}
	}
}

}

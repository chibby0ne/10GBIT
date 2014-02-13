//
//  Copyright (C) 2010 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Hardware-compliant LDPC decoder shared functionality
/// \author Matthias Alles
/// \date   2010/12/03
//

#include "dec_ldpc_bin_hw_share.h"

/*
 * 0 : no logging
 * 1 : results on iteration level (parity checks failed)
 * 2 : results on check node group level (parity checks failed)
 * 3 : check node inputs / outputs
 */
#define LOG_LEVEL_LDPC 3
#define NLOGGING
#include <iostream>
#include "../../cse/assistance/debug_macros.h"
#include <cassert>


using namespace hlp_fct::math;
using std::cout;
using std::endl;

namespace cse_lib {

template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Change_Shift_Direction()
{
	for (unsigned int i = 0; i < max_check_degree_ * num_check_nodes_ / dst_parallelism_; i++)
		if (shft_vector_[i] > 0)
			shft_vector_[i] = dst_parallelism_ - shft_vector_[i];
}


template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Init_APP_RAM(bool          parity_reordering,
                                                   Buffer<T>    &input_bits_llr,
                                                   Buffer<T, 2> &app_ram)
{
	unsigned int folding_factor;
	unsigned int num_info_nodes;
	unsigned int max_addr_normal_distribution;
	unsigned int i, j;

	folding_factor = src_parallelism_ / dst_parallelism_;
	num_info_nodes = num_variable_nodes_ - num_check_nodes_;

	if (parity_reordering) {

		// For IRA codes we only distribute the information nodes vertically
		max_addr_normal_distribution = num_info_nodes / dst_parallelism_;

	} else {

		// For structured codes all bits are distributed the same
		max_addr_normal_distribution = num_variable_nodes_ / dst_parallelism_;

	}

	/*
	 * For each address distribute the channel values according to the decoder
	 * parallelism. For IRA codes this is done only for the systematic part
	 * like that.
	 */
	for(i = 0; i < max_addr_normal_distribution; i++)
		for(j = 0; j < dst_parallelism_; j++)
			app_ram[j][i] = input_bits_llr[
			                ((i / folding_factor) * dst_parallelism_ * folding_factor) +
			                 (j * folding_factor) +  // folding factor steps wide
			                 (i % folding_factor)];  // offset within submatrix

	if (parity_reordering) {
	
		// The parity nodes are distributed vertically for IRA codes!
		for(j = 0 ; j < dst_parallelism_; j++)
			for(i = 0; i < num_check_nodes_ / dst_parallelism_; i++)
				app_ram[j][i + num_info_nodes / dst_parallelism_] =
					    input_bits_llr[num_info_nodes +
					              i / folding_factor +
					              (i % folding_factor) * num_check_nodes_ / src_parallelism_ +
					              j * num_check_nodes_ / src_parallelism_ * folding_factor];
	}

	return;
}


template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Read_APP_RAM(Buffer<T, 2>            &app_ram,
                                                   int                      iter,
                                                   Buffer<T, 2>            &output_bits_llr_app,
                                                   Buffer<unsigned int, 2> &output_bits)
{
	unsigned int folding_factor;
	unsigned int num_info_nodes;
	unsigned int max_addr_normal_distribution;
	unsigned int i, j;

	folding_factor = src_parallelism_ / dst_parallelism_;
	num_info_nodes = num_variable_nodes_ - num_check_nodes_;


	if (is_IRA_code_) {

		// For IRA codes we only read the information nodes vertically
		max_addr_normal_distribution = num_info_nodes / dst_parallelism_;

	} else {

		// For structured codes all bits are distributed the same
		max_addr_normal_distribution = num_variable_nodes_ / dst_parallelism_;

	}

	/*
	 * For each address read the channel values according to the decoder
	 * parallelism. For IRA codes this is done only for the systematic part
	 * like that.
	 */
	for(i = 0; i < max_addr_normal_distribution; i++)
		for(j = 0; j < dst_parallelism_; j++)
			output_bits_llr_app[iter][ ((i / folding_factor) * dst_parallelism_ * folding_factor) +
			                                     (j * folding_factor) +
			                                     (i % folding_factor)] = app_ram[j][i];

	if (is_IRA_code_) {

		// The parity nodes are distributed vertically for IRA codes!
		for(j = 0 ; j < dst_parallelism_; j++)
			for(i = 0; i < num_check_nodes_ / dst_parallelism_; i++)
				output_bits_llr_app[iter][num_info_nodes +
				            i / folding_factor +
				           (i % folding_factor) * num_check_nodes_ / src_parallelism_ +
				            j * num_check_nodes_ / src_parallelism_ * folding_factor] =
				          app_ram[j][i + num_info_nodes / dst_parallelism_];
	}

	// Calculate the sign bit.
	for(i = 0; i < num_variable_nodes_ ; i++)
		output_bits[iter][i] = (output_bits_llr_app[iter][i] < 0 ? 1 : 0);

	return;

}


template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Get_Check_Node_Input(Buffer<T, 2> &app_ram,
                                                           Buffer<T, 2> &msg_ram,
                                                           unsigned int  iter,
                                                           unsigned int  cng_counter,
                                                           unsigned int  cfu_counter,
                                                           Buffer<T>    &check_node_in,
                                                           Buffer<T>    &app_out)
{
	int vng_select;
	int shift_value;
	T app_ram_content;
	unsigned int vn_select;
	unsigned int vector_addr;

	// Calculate the base address for addressing shift and address vectors.
	vector_addr = cng_counter * max_check_degree_;

	// Iterate over the single edges.
	for(unsigned int vn2cn_msg = 0; vn2cn_msg < max_check_degree_; vn2cn_msg++)
	{
		// Determine the variable node group we want to process and its shift value.
		vng_select = addr_vector_[vector_addr];
		shift_value = shft_vector_[vector_addr];

		// Determine the variable node within the current variable node group.
		vn_select = (shift_value + cfu_counter) % dst_parallelism_;

		/*
		 * Calulate the input for the check nodes here. Only if the
		 * address to the variable node group is valid, we have to
		 * read the app_ram. Otherwise we insert the neutral element
		 * with respect to the check node operation (= maximum
		 * message).
		 */
		if(vng_select > -1)
		{
			app_ram_content = app_ram[vn_select][vng_select];

			// Calculate check node input.
			if(iter != 0)

				// Subtract extrinsic information from previous iteration.
				check_node_in[vn2cn_msg] = app_ram_content - msg_ram[vn_select][vector_addr];

			else

				// In the first iteration read only the channel information.
				check_node_in[vn2cn_msg] = app_ram_content;

			app_out[vn2cn_msg] = app_ram_content;
		}
		else
		{
			// No valid address => Insert value that does not affect the check node.
			//check_node_in[vn2cn_msg] = numeric_limits<T>::max();
			check_node_in[vn2cn_msg] = max_msg_extr_;

			app_out[vn2cn_msg] = 0;
		}

		if(is_IRA_code_     &&
		   cfu_counter == 0 &&
		   shift_value == (signed) dst_parallelism_ - 1 &&
		   vng_select  == (signed) (num_variable_nodes_ / dst_parallelism_ - 1))
		{
			//check_node_in[vn2cn_msg] = numeric_limits<T>::max();
			check_node_in[vn2cn_msg] = static_cast<int> (max_msg_extr_);
			app_out[vn2cn_msg] = 0;
		}

		/*
		 * Saturate the input to the check node, symmetrically due to signed-magnitude
		 * representation in hardware. For float models, we do not saturate here.
		 */
		Saturate(check_node_in[vn2cn_msg], max_msg_extr_);

		// Update to next edge.
		vector_addr++;
	}

	LM_OUT_LEVEL(LDPC, 3, "CNG: " << cng_counter <<
	                      " CFU: " << cfu_counter <<
	                      " APP RAM shifted out: " << app_out << endl);

}


template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Write_Check_Node_Output(Buffer<T, 2> &app_ram,
                                                              Buffer<T, 2> &msg_ram,
                                                              unsigned int  iter,
                                                              unsigned int  cng_counter,
                                                              unsigned int  cfu_counter,
                                                              Buffer<T>    &check_node_out)
{
	int vng_select;
	int shift_value;
	unsigned int vn_select;
	unsigned int vector_addr;

	// Calculate the base address for shift and address vectors.
	vector_addr = cng_counter * max_check_degree_;

	// Process each edge of the check node separately.
	for(unsigned int cn2vn_msg = 0; cn2vn_msg < max_check_degree_; cn2vn_msg++)
	{
		// Determine the variable node group we want to process.
		vng_select = addr_vector_[vector_addr];
		shift_value = shft_vector_[vector_addr];

		// Determine the variable node within the current variable node group.
		vn_select = (shift_value + cfu_counter) % dst_parallelism_;

		// In case of IRA codes set the output for virtual edge to 0.
		if(is_IRA_code_     &&
		   cfu_counter == 0 &&
		   shift_value == (signed) dst_parallelism_ - 1 &&
		   vng_select  == (signed) (num_variable_nodes_ / dst_parallelism_ - 1))
		{
			check_node_out[cn2vn_msg] = 0;
		}

		// In case of a valid variable node group, update app_ram and msg_ram.
		if(vng_select > -1)
		{
			T current_app_value = app_ram[vn_select][vng_select];
			T current_message   = check_node_out[cn2vn_msg];

			if (iter != 0)
				current_app_value -= msg_ram[vn_select][vector_addr];

			// Saturate already here in order to be 100% hardware-compliant!
			Saturate(current_app_value, max_msg_app_, - max_msg_app_ - 1);

			current_app_value += current_message;

			/*
			 * Saturate the APP value. In contrast to the input check node messages,
			 * there is no need to saturate symmetrically, since this value will
			 * never be represented as a signed-magnitude value.
			 */
			Saturate(current_app_value, max_msg_app_, - max_msg_app_ - 1);

			// Write back results
			app_ram[vn_select][vng_select]  = current_app_value;
			msg_ram[vn_select][vector_addr] = current_message;
		}

		vector_addr++;
	}
}


template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Calc_Parity_Check(Buffer<T> &values)
{

	unsigned int parity_check = 0;

	for(unsigned int i = 0; i < max_check_degree_; i++)
		if (values[i] < 0)
			parity_check++;

	parity_check %= 2;

	return parity_check;
}


/*********************
 * Decoder Functions *
 *********************/


template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Decode_Layered(Buffer<T, 2> &app_ram,
                                                             Buffer<T, 2> &msg_ram,
                                                             int           iter)
{
	Buffer<T> check_node_io(max_check_degree_);
	Buffer<T> app_out(max_check_degree_);
	cn_msg_abs_.Resize(max_check_degree_);
	cn_msg_sign_.Resize(max_check_degree_);
	unsigned int ok_checks = 0;

	// Iterate over all check node groups.
	for(unsigned int cng_counter = 0; cng_counter < num_check_nodes_ / dst_parallelism_; cng_counter++)
	{
		// Iterate over the functional units.
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
		{
			Get_Check_Node_Input(app_ram, msg_ram, iter, cng_counter, cfu_counter, check_node_io, app_out);

			ok_checks += (1 - Calc_Parity_Check(app_out));

			switch(check_node_algorithm_)
			{
            case Decoder_LDPC_Check_Node_Share::MIN_SUM_SELF_CORRECTING:
			case Decoder_LDPC_Check_Node_Share::MIN_SUM:
				Check_Node_Min_Sum(check_node_io, esf_factor_); break;

			case Decoder_LDPC_Check_Node_Share::LAMBDA_MIN:
				Check_Node_Lambda_Min(check_node_io, num_lambda_min_, bw_fract_); break;

			case Decoder_LDPC_Check_Node_Share::LAMBDA_3MIN_3MAG:
				Check_Node_Lambda_3Min_3Mag(check_node_io, bw_fract_); break;

            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW:
            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW_IMPROVED:
            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW_SELF_CORRECTING:
                Check_Node_Split_Row(check_node_io, esf_factor_, num_partitions_, threshold_); break;

			}

			Write_Check_Node_Output(app_ram, msg_ram, iter, cng_counter, cfu_counter, check_node_io);
		}
	}

	return ok_checks;
}


template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Decode_Layered_Superposed(Buffer<T, 2> &app_ram,
                                                                        Buffer<T, 2> &msg_ram,
                                                                        int          iter)
{
	Buffer<T, 2> check_node_io(dst_parallelism_, max_check_degree_);
	Buffer<T>    app_out(max_check_degree_);
	cn_msg_abs_.Resize(max_check_degree_);
	cn_msg_sign_.Resize(max_check_degree_);
	unsigned int ok_checks = 0;

	// Iterate over all check node groups.
	for(unsigned int cng_loop_counter = 0; cng_loop_counter < num_check_nodes_ / dst_parallelism_; cng_loop_counter++)
	{
		unsigned int cng_counter;

		//if (iter % 2 == 0)
			cng_counter = cng_loop_counter;
				//else
			//cng_counter = num_check_nodes_ / dst_parallelism_ - 1 - cng_loop_counter;

		// Read the input for a whole check node group first.
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
		{
			Get_Check_Node_Input(app_ram,
			                     msg_ram,
			                     iter,
			                     cng_counter,
			                     cfu_counter,
			                     check_node_io[cfu_counter],
			                     app_out);

			ok_checks += (1 - Calc_Parity_Check(app_out));
		}

		LM_OUT_LEVEL(LDPC, 2, "CNG: " << cng_counter <<
		                      " Checks failed: " << dst_parallelism_ * (1 + cng_counter) - ok_checks << endl);

		// Perform the Check Node calculations for all CFUs now.
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
		{
			LM_OUT_LEVEL(LDPC, 3, "CNG: " << cng_counter <<
			                      " CFU: " << cfu_counter <<
			                      " in: " << check_node_io[cfu_counter] << endl);

			switch(check_node_algorithm_)
			{
			case Decoder_LDPC_Check_Node_Share::MIN_SUM:
            case Decoder_LDPC_Check_Node_Share::MIN_SUM_SELF_CORRECTING:
				Check_Node_Min_Sum(check_node_io[cfu_counter], esf_factor_);
				break;

			case Decoder_LDPC_Check_Node_Share::LAMBDA_MIN:
				Check_Node_Lambda_Min(check_node_io[cfu_counter], num_lambda_min_, bw_fract_);
				break;

			case Decoder_LDPC_Check_Node_Share::LAMBDA_3MIN_3MAG:
				Check_Node_Lambda_3Min_3Mag(check_node_io[cfu_counter], bw_fract_);
				break;

            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW:
            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW_IMPROVED:
            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW_SELF_CORRECTING:
                Check_Node_Split_Row(check_node_io[cfu_counter], esf_factor_, num_partitions_, threshold_);
                break;
			}

			LM_OUT_LEVEL(LDPC, 3, "CNG: " << cng_counter <<
			                      " CFU: " << cfu_counter <<
			                      " out: " << check_node_io[cfu_counter] << endl);
		}

		// Write back the result of all check nodes.
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
			Write_Check_Node_Output(app_ram,
			                        msg_ram,
			                        iter,
			                        cng_counter,
			                        cfu_counter,
			                        check_node_io[cfu_counter]);
	}

	LM_OUT_LEVEL(LDPC, 1, "Iter: " << iter << " Checks failed: " << num_check_nodes_ - ok_checks << endl);

	return ok_checks;
}


template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Decode_Two_Phase(Buffer<T, 2> &app_ram,
                                                               Buffer<T, 2> &msg_ram,
                                                               int           iter,
                                                               bool          app_parity_check)
{
	// Number of check node groups
	unsigned int num_cng = num_check_nodes_ / dst_parallelism_;

	Buffer<T, 3> check_node_io;
	unsigned int dimensions[3];
	dimensions[0] = num_cng;
	dimensions[1] = dst_parallelism_;
	dimensions[2] = max_check_degree_;
	check_node_io.Resize(dimensions);
	Buffer<T> app_out(max_check_degree_);
	cn_msg_abs_.Resize(max_check_degree_);
	cn_msg_sign_.Resize(max_check_degree_);
	unsigned int ok_checks = 0;

	// Iterate over all check nodes.
	for(unsigned int cng_counter = 0; cng_counter < num_cng; cng_counter++)
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
		{
			Get_Check_Node_Input(app_ram,
			                     msg_ram,
			                     iter,
			                     cng_counter,
			                     cfu_counter,
			                     check_node_io[cng_counter][cfu_counter],
			                     app_out);

			if (app_parity_check)
				ok_checks += (1 - Calc_Parity_Check(app_out));
			else
				ok_checks += (1 - Calc_Parity_Check(check_node_io[cng_counter][cfu_counter]));
		}

	/*
	 * Break current iteration if the all parity checks are successful.
	 * Since the current iteration may change the decision,
	 * do not write the app values back
	 */
	if (ok_checks == num_cng * dst_parallelism_ && app_parity_check)
		return ok_checks;

	// Perform the Check Node calculations for all CFUs now.
	for(unsigned int cng_counter = 0; cng_counter < num_cng; cng_counter++)
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
			switch(check_node_algorithm_)
			{
			case Decoder_LDPC_Check_Node_Share::MIN_SUM:
            case Decoder_LDPC_Check_Node_Share::MIN_SUM_SELF_CORRECTING:
				Check_Node_Min_Sum(check_node_io[cng_counter][cfu_counter], esf_factor_);
				break;

			case Decoder_LDPC_Check_Node_Share::LAMBDA_MIN:
				Check_Node_Lambda_Min(check_node_io[cng_counter][cfu_counter], num_lambda_min_, bw_fract_);
				break;

			case Decoder_LDPC_Check_Node_Share::LAMBDA_3MIN_3MAG:
				Check_Node_Lambda_3Min_3Mag(check_node_io[cng_counter][cfu_counter], bw_fract_);
				break;
			case Decoder_LDPC_Check_Node_Share::SPLIT_ROW:
            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW_IMPROVED:
            case Decoder_LDPC_Check_Node_Share::SPLIT_ROW_SELF_CORRECTING:
				Check_Node_Split_Row(check_node_io[cng_counter][cfu_counter], esf_factor_, num_partitions_, threshold_);
				break;
			}

	// Write back the result of all check nodes.
	for(unsigned int cng_counter = 0; cng_counter < num_cng; cng_counter++)
		for(unsigned int cfu_counter = 0; cfu_counter < dst_parallelism_; cfu_counter++)
			Write_Check_Node_Output(app_ram,
			                        msg_ram,
			                        iter,
			                        cng_counter,
			                        cfu_counter,
			                        check_node_io[cng_counter][cfu_counter]);

	return ok_checks;
}



/************************
 * Check Node Functions *
 ************************/

template<class T>
int Decoder_LDPC_Binary_HW_Share<T>::Check_Node_Lambda_Min(Buffer<T>    &in_out_msg,
                                                        unsigned int  num_lambda_min,
                                                        unsigned int  num_bits_fract)
{
	unsigned int check_node_degree = in_out_msg.length();
	unsigned int min_msg[4];
	unsigned int min_idx[4];
	int pchk_sign = 1;
	unsigned int minstar_msg;
	unsigned int index = 0;

	// Separate magnitude and sign.
	for (unsigned int i = 0; i < check_node_degree; i++)
	{
		cn_msg_abs_[i]  = abs(in_out_msg[i]);
		cn_msg_sign_[i] = (in_out_msg[i] >= 0) ? 1 : -1;
		pchk_sign *= cn_msg_sign_[i];
	}

	// Find the num_lambda_min() minima.
	for(unsigned int i = 0; i < num_lambda_min; i++)
	{
		// Find the current minimum and remember its magnitude and index.
		min_msg[i] = Minimum_Fixp(cn_msg_abs_, index);
		min_idx[i] = index;

		// Replace the current minimum and continue with next minimum search.
//		cn_msg_abs_[index] = numeric_limits<T>::max();
		cn_msg_abs_[index] = max_msg_extr_ + 1;
	}

	// 2-Min algorithm
	if (num_lambda_min == 2)
		for(unsigned int i = 0; i < check_node_degree; i++)
		{
			if(min_idx[0] == i)
				minstar_msg = min_msg[1];

			else if(min_idx[1] == i)
				minstar_msg = min_msg[0];

			else
				minstar_msg = Minstar_Fixp(min_msg[0], min_msg[1], num_bits_fract);

			// Calculate new outgoing edge with new sign and new magnitude.
			in_out_msg[i] = minstar_msg * cn_msg_sign_[i] * pchk_sign;
		}

	// 3-Min algorithm
	else if (num_lambda_min == 3)
		for(unsigned int i = 0; i < check_node_degree; i++)
		{
			if(min_idx[0] == i)
				minstar_msg = Minstar_Fixp(min_msg[1], min_msg[2], num_bits_fract);

			else if(min_idx[1] == i)
				minstar_msg = Minstar_Fixp(min_msg[0], min_msg[2], num_bits_fract);

			else if(min_idx[2] == i)
				minstar_msg = Minstar_Fixp(min_msg[0], min_msg[1], num_bits_fract);

			else
				minstar_msg = Minstar_Fixp(Minstar_Fixp(min_msg[0], min_msg[1], num_bits_fract),
				                           min_msg[2],
				                           num_bits_fract);

			// Calculate new outgoing edge with new sign and new magnitude.
			in_out_msg[i] = minstar_msg * cn_msg_sign_[i] * pchk_sign;
		}

	// Return whether parity check was satisfied or not.
	return pchk_sign;
}


template<class T>
int Decoder_LDPC_Binary_HW_Share<T>::Check_Node_Min_Sum(Buffer<T> &in_out_msg, float esf_factor)
{
	unsigned int check_node_degree = in_out_msg.length();
	T min0_msg, min1_msg, min_out_msg;
	unsigned int min0_idx;
	int  pchk_sign = 1;
	unsigned int index = 0;

	// Separate magnitude and sign.
	for (unsigned int i = 0; i < check_node_degree; i++)
	{
		cn_msg_abs_[i]  = abs(in_out_msg[i]);
		cn_msg_sign_[i] = (in_out_msg[i] >= 0) ? 1 : -1;
		pchk_sign *= cn_msg_sign_[i];
	}

	// Find the current minimum and remember its magnitude and index.
	min0_msg = Minimum_Fixp(cn_msg_abs_, index);
	min0_idx = index;

	// Replace the current minimum and continue with next minimum search.
//	cn_msg_abs_[index] = numeric_limits<T>::max();
	cn_msg_abs_[index] = max_msg_extr_;

	// Find the second minimum and remember its magnitude and index.
	min1_msg = Minimum_Fixp(cn_msg_abs_, index);

	Scale_ESF(min0_msg, esf_factor);
	Scale_ESF(min1_msg, esf_factor);

	for(unsigned int i = 0; i < check_node_degree; i++)
	{
		if(min0_idx == i)
			min_out_msg = min1_msg;
		else
			min_out_msg = min0_msg;

		// Calculate new outgoing edge with new sign and new magnitude.
		in_out_msg[i] = min_out_msg * cn_msg_sign_[i] * pchk_sign;
	}

	return pchk_sign;
}

// Check_Node_Split_Row: implements Check_Node functionality according to Split-Row Treshold and Split-Row Threshold Improved
template<class T>
int Decoder_LDPC_Binary_HW_Share<T>::Check_Node_Split_Row(Buffer<T> &in_out_msg, float esf_factor, unsigned int partitions, int threshold)
{

    std::vector< std::vector<T> > parts;      // vector that holds partitions
    std::vector< std::vector<T> > temp_parts;      // vector that holds partitions temporarily
    parts.resize(partitions);
    temp_parts.resize(partitions);
    int sign_part[partitions];               // array that holds accumulated sign of each partition
    T local_min1[partitions];     // array that holds min1 for each partition
    T local_min2[partitions];     // array that holds min2 for each partition
    unsigned int min1_idx[partitions];       // array that hold index of min1 for each partition
    std::vector<unsigned int> threshold_en(partitions, 0);   // vector that holds the threshold_en signal for each partition
    unsigned int len = in_out_msg.length() / partitions;     // length of each partition = number of edges per partition
    int pchk_sign = 1;


   // distribute the messages between the partitions
   Distribute_Messages_Partitions(parts, in_out_msg, partitions, len);



   // find local first and second minimum, signs and set threshold_enable for each partition
   for (unsigned int part_num = 0; part_num != partitions; part_num++) {
   sign_part[part_num] = Split_Row_Local_Minimum_Sign(parts[part_num], local_min1[part_num], min1_idx[part_num], local_min2[part_num], threshold_en[part_num], threshold, part_num);
   }


   // apply Split Row Threshold Improved or Split Row Threshold (Split Row Self Correcting as well) algorithm
   if (check_node_algorithm_ == SPLIT_ROW_IMPROVED || check_node_algorithm_ == SPLIT_ROW_SELF_CORRECTING) {

       for (unsigned int part_num = 0; part_num != partitions; part_num++) {
           if (local_min1[part_num] <= threshold && local_min2[part_num] <= threshold) {
               // condition 1
               Update_Minimum(temp_parts[part_num], min1_idx[part_num], local_min1[part_num], local_min2[part_num], len);
           } else if (local_min1[part_num] <= threshold && local_min2[part_num] > threshold) {
               // condition 2
               if (Neighbors_Threshold_Enable(part_num, partitions, threshold_en)) {
                   // condition 2a
                   Update_Minimum(temp_parts[part_num], min1_idx[part_num], local_min1[part_num], threshold, len);
               } else {
                   // condition 2b
                   Update_Minimum(temp_parts[part_num], min1_idx[part_num], local_min1[part_num], local_min2[part_num], len);
               }
           } else if (local_min1[part_num] > threshold && Neighbors_Threshold_Enable(part_num, partitions, threshold_en)) {
               // condition 3
               Update_Minimum(temp_parts[part_num], min1_idx[part_num], threshold, threshold, len);
           } else {
               // condition 4
               Update_Minimum(temp_parts[part_num], min1_idx[part_num], local_min1[part_num], local_min2[part_num], len);
           }
           Multiply_Sign_Esf(parts[part_num], temp_parts[part_num], sign_part[part_num], esf_factor, part_num);
       }

   } else {

       for (unsigned int part_num = 0; part_num != partitions; part_num++) {
           if (local_min1[part_num] > threshold && Neighbors_Threshold_Enable(part_num, partitions, threshold_en)) {
               // condition 3
               Update_Minimum(temp_parts[part_num], min1_idx[part_num], threshold, threshold, len);
           } else {
               // condition 1 and 4
               Update_Minimum(temp_parts[part_num], min1_idx[part_num], local_min1[part_num], local_min2[part_num], len);
           }
           // multiply min with local sign and own sign
           Multiply_Sign_Esf(parts[part_num], temp_parts[part_num], sign_part[part_num], esf_factor, part_num);
       }
   }

   // Also consider the sign of the adjacent partitions
    for (unsigned int part_num = 1; part_num != partitions; part_num++) {
        Multiply_Split_Sign_Signal(parts[part_num - 1], sign_part[part_num]);
        Multiply_Split_Sign_Signal(parts[part_num], sign_part[part_num - 1]);
        if (part_num == 1) {
            // if is the first time, then multiply the sign of first and second partition
            // else keep multiply new partition sign with the accumulated one
            pchk_sign = sign_part[part_num - 1] * sign_part[part_num];
        } else {
            pchk_sign *= sign_part[part_num];
        }
    }

    // put everything in in_out_msg so it can be read in the variable node again
    Regroup_Messages_Partitions(parts, in_out_msg, len);

    return pchk_sign;
}

// Distribute_Messages_Partitions: distributes messages from variable node (in_out_msg) between partitions
template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Distribute_Messages_Partitions(std::vector< std::vector<T> > &parts, Buffer<T> &in_out_msg, unsigned int partitions, unsigned int len)
{
    for (unsigned int i = 0; i != partitions; i++) {
        for (unsigned int j = 0; j != len; j++) {
            parts[i].push_back(in_out_msg[i * len + j]);
        }
    }
}


//int Decoder_LDPC_Binary_HW_Share::Split_Row_Local_Minimum_Sign(std::vector<int>& part, unsigned int& minimum1, unsigned int& min1_idx, unsigned int& minimum2, unsigned int& threshold_en, unsigned int threshold, unsigned int part_num)
// Split_Row_Local_Minimum_Sign: gets local minima (1st and 2nd), gets index of 1st minima, sets threshold_en signal, and calculates split_row_sign signal of all msg per partition
template<class T>
int Decoder_LDPC_Binary_HW_Share<T>::Split_Row_Local_Minimum_Sign(std::vector<T>& part, T& minimum1, unsigned int& min1_idx, T& minimum2, unsigned int& threshold_en, int threshold, unsigned int part_num)
{
    unsigned int check_node_degree = part.size();
    assert(check_node_degree == 16);
	int  pchk_sign = 1;
	unsigned int index = 0;

	// Separate magnitude and sign.
	for (unsigned int i = 0; i != check_node_degree; i++)
	{
        // we need part_num in order to store them in a different place from other partitions (remember multiple partitions)
        // we are storing all magnitudes and signs of each edge in one buffer, but we use indexing to access correct one for each partition
        // part does not need additional part_num * check_node_degree because this is the vector for that partition (vector<vector<unsigned int> > is taken care of in caller) 
		cn_msg_abs_[part_num * check_node_degree + i]  = abs(part[i]);
		cn_msg_sign_[part_num * check_node_degree + i] = (part[i] >= 0) ? 1 : -1;
		pchk_sign *= cn_msg_sign_[part_num * check_node_degree + i];
	}

	// Find the current minimum and remember its magnitude and index.
	minimum1 = Minimum_Fixp_Split_Row(cn_msg_abs_, index, part_num, check_node_degree);
	min1_idx = index; // Update_Min searches in that partition from 0 to len = (in_out_msg.length() / partitions)

	// Replace the current minimum and continue with next minimum search.
	cn_msg_abs_[part_num * check_node_degree + index] = max_msg_extr_;

	// Find the second minimum and remember its magnitude and index.
	minimum2 = Minimum_Fixp_Split_Row(cn_msg_abs_, index, part_num, check_node_degree);

    // store minimums
    // in the case of split row, this is done in function Update_minimum

    // if local min <= threshold, then threshold_en = 1
    if (minimum1 <= threshold) {
        threshold_en = true;
    } else {
        threshold_en = false;
    }
    return pchk_sign;
}


//bool Decoder_LDPC_Binary_HW_Share::Neighbors_Threshold_Enable(unsigned int part_num, unsigned int partitions, unsigned int threshold_en[])
// Neighbors_Threshold_Enable: returns true if any of the neigboring partitions has threshold_en asserted, else false
template<class T>
bool Decoder_LDPC_Binary_HW_Share<T>::Neighbors_Threshold_Enable(unsigned int part_num, unsigned int partitions, std::vector<unsigned int>& threshold_en)
{
    if (part_num == 0) {
        return threshold_en[1];
    } else if (part_num == partitions - 1) {
        return threshold_en[part_num - 1];
    } else {
        return threshold_en[part_num - 1] || threshold_en[part_num + 1];
    }
}

// Update_Minimum: stores the magnitude of that splitted row
template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Update_Minimum(std::vector<T>& temp_part, unsigned int min1_idx, T min1, T min2, unsigned int len)
{
    T min_out_msg;
    unsigned int check_node_degree = len;

    for (unsigned int i = 0; i != check_node_degree; i++) {
        if (i == min1_idx) {
            min_out_msg = min2;
        } else {
            min_out_msg = min1;
        }
        temp_part.push_back(min_out_msg);
    }
}

// Multiply_Sign_Esf: Multiplies magnitude and local sign and esf factor and stores it in signed representation (vector <> part)
template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Multiply_Sign_Esf(std::vector<T>& part, std::vector<T>& temp_part, int local_pchk_sign, float esf_factor, unsigned int part_num)
{
    assert(part.size() == temp_part.size());
    unsigned int check_node_degree = part.size();
	for(unsigned int i = 0; i < check_node_degree; i++)
	{
        Scale_ESF(temp_part[i], esf_factor);

        // Calculate new outgoing edge with new sign and new magnitude.
        // we need j in order to multiply it with correct sign (remember multiple partitions)
		part[i] = temp_part[i] * cn_msg_sign_[part_num * check_node_degree + i] * local_pchk_sign;
	}
}

// Multiply_Split_Sign_Signal: Multiply messages of each partition with the sign of the neighboring partition and store it as message
template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Multiply_Split_Sign_Signal(std::vector<T>& part, int split_sign)
{
    unsigned int check_node_degree = part.size();

    for (unsigned int i = 0; i != check_node_degree; i++) {
        part[i] *= split_sign;
    }
}

// Regroup_Messages_Partitions: Regroups messages of all partitions in Buffer<in_out_msg> used in Write_Check_Node_Output
template<class T>
void Decoder_LDPC_Binary_HW_Share<T>::Regroup_Messages_Partitions(std::vector< std::vector<T> > &parts, Buffer<T> &in_out_msg, unsigned int len)
{
    for (unsigned int i = 0, j = 0, k = 0; i != in_out_msg.length(); i++, k++) {
        /* 
         * if j reaches the lenght of a row of a partition
         * then we have finished reading this partition,
         * so we grab the next partition
         */
        if (i != 0 && i % len == 0) {
            j++;
            k = 0;
        }
        in_out_msg[i] = parts[j][k];
    }

}

// Minimum_Fixp_Split_Row: finds minima and its index in the messages of a partition
template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Minimum_Fixp_Split_Row(Buffer<T> &message, unsigned int& idx, unsigned int part_num, unsigned int check_node_degree)
{
    T temp = max_msg_extr_;

    // here check_node_degree = len = length of a partition
    for(unsigned int i = 0; i != check_node_degree; i++) {
        if (message[part_num * check_node_degree + i] <= temp) {
            idx  = i;
            temp = message[part_num * check_node_degree + i];
        }
    }

    return temp;
}


template<class T>
int Decoder_LDPC_Binary_HW_Share<T>::Check_Node_Lambda_3Min_3Mag(Buffer<T> &in_out_msg,
                                                              unsigned int num_bits_fract)
{
	unsigned int check_node_degree = in_out_msg.length();
	unsigned int min_msg[4];
	unsigned int min_idx[4];
	int pchk_sign = 1;
	unsigned int minstar_msg;
	unsigned int index = 0;

	// Separate magnitude and sign.
	for (unsigned int i = 0; i < check_node_degree; i++)
	{
		cn_msg_abs_[i]  = abs(in_out_msg[i]);
		cn_msg_sign_[i] = (in_out_msg[i] >= 0) ? 1 : -1;
		pchk_sign *= cn_msg_sign_[i];
	}

	// Find the 3 minima.
	for(unsigned int i = 0; i < 3; i++)
	{
		// Find the current minimum and remember its magnitude and index.
		min_msg[i] = Minimum_Fixp(cn_msg_abs_, index);
		min_idx[i] = index;

		// Replace the current minimum and continue with next minimum search.
//		cn_msg_abs_[index] = numeric_limits<T>::max();
		cn_msg_abs_[index] = max_msg_extr_;
	}

	unsigned int minstar_msg_else_branch = Minstar_Fixp(min_msg[0], min_msg[1], num_bits_fract);

	for (unsigned int i = 0; i < check_node_degree; i++)
	{
		if (min_idx[0] == i)
			minstar_msg = Minstar_Fixp(min_msg[1], min_msg[2], num_bits_fract);

		else if (min_idx[1] == i)
			minstar_msg = Minstar_Fixp(min_msg[0], min_msg[2], num_bits_fract);

		else
			minstar_msg = minstar_msg_else_branch;

		// Calculate new outgoing edge with new sign and new magnitude.
		in_out_msg[i] = minstar_msg * cn_msg_sign_[i] * pchk_sign;
	}

	// Return whether parity check was satisfied or not.
	return pchk_sign;
}


template<class T>
T Decoder_LDPC_Binary_HW_Share<T>::Minimum_Fixp(Buffer<T>    &message,
                                                unsigned int &idx)
{
	T temp = message[0];
	idx    = 0;

	for(unsigned int i = 1; i < message.length(); i++)
		if(message[i] <= temp)
		{
			idx  = i;
			temp = message[i];
		}

	return temp;
}


template<class T>
inline int Decoder_LDPC_Binary_HW_Share<T>::Delta(int x, unsigned int num_bits_fract)
{

	// Hardware-compliant function look-up.
	if (num_bits_fract == 2)
	{
		if (x < 8)
			if (x < 4)
				return 2;
			else
				return 1;
		else
			return 0;
	}
	else
	{
		// generic function
		int d = ((5 << (num_bits_fract + 2 - 3)) - (x)) >> 2;

		return (d > 0) ? d : 0;
	}
}


template<class T>
inline float Decoder_LDPC_Binary_HW_Share<T>::Delta(float x, unsigned int)
{
	return log ( 1 + exp ( - x));
}


template<class T>
inline int Decoder_LDPC_Binary_HW_Share<T>::Minstar_Fixp(T a, T b, unsigned int num_bits_fract)
{
	int min = (a < b) ? a : b;
	int d = - Delta( abs(a - b), num_bits_fract ) + Delta(a + b, num_bits_fract);

	return (min + d);
}


template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Calc_Modified_Systematic_Bits(unsigned int             iter,
                                                                            Buffer<T>               &input_bits_llr,
                                                                            Buffer<unsigned int, 2> &output_bits)
{
	unsigned int num_modified_systematic_bits = 0;

	// Calculate the number of modified signs of the systematic information.
	for(unsigned int i = 0; i < num_variable_nodes_ - num_check_nodes_; i++)
	{
		/*
		 * If input was 1 (negative LLR) and output is 0 or
		 * if input was 0 (positive LLR) and output is 1 then
		 * increment number of modified systematic bits.
		 */
		if((output_bits[iter - 1][i] == 0 && input_bits_llr[i]  < 0) ||
		   (output_bits[iter - 1][i] == 1 && input_bits_llr[i] >= 0))
		{
			num_modified_systematic_bits++;
		}
	}

	return num_modified_systematic_bits;
}


template<class T>
inline void Decoder_LDPC_Binary_HW_Share<T>::Scale_ESF(int &msg, float esf_factor)
{
	if(esf_factor == 0.875)
	{
		// Multiply with 0.875 like it is done in the hardware.
		msg = ( ((msg << 1) +
		         (msg << 0) +
		         (msg >> 1) + 1) >> 2);
	}
	else
	{
		// Multiply with 0.75 like it is done in the hardware.
		msg = ( ((msg << 1) + (msg << 0) + 1) >> 2);
	}
}


template<class T>
inline void Decoder_LDPC_Binary_HW_Share<T>::Scale_ESF(float &msg, float esf_factor)
{
	msg *= esf_factor;
}


template<class T>
inline void Decoder_LDPC_Binary_HW_Share<T>::Saturate(int& input, int max_msg)
{
	Saturate_Value(input, max_msg);
}


template<class T>
inline void Decoder_LDPC_Binary_HW_Share<T>::Saturate(float&, float)
{
}


template<class T>
inline void Decoder_LDPC_Binary_HW_Share<T>::Saturate(int &input, int min_msg, int max_msg)
{
	Saturate_Value(input, min_msg, max_msg);
}


template<class T>
inline void Decoder_LDPC_Binary_HW_Share<T>::Saturate(float&, float, float)
{
}

template<class T>
unsigned int Decoder_LDPC_Binary_HW_Share<T>::Calc_Flipped_Bits(unsigned int iter, Buffer<unsigned int, 2> &output_bits)
{
    unsigned int flipped_bits = 0;
    // Calculate the number of flipped (sign) bits between iterations
    for (unsigned int i = 0; i != num_variable_nodes_; i++) {
        /*
         * If sign bit in current iteration is different from sign bit from previous
         * then increment number of flipped bits
         */
        if (output_bits[iter - 1][i] != output_bits[iter][i]) {
            flipped_bits++;
        }
    }
    return flipped_bits;
}


template class Decoder_LDPC_Binary_HW_Share<int>;
template class Decoder_LDPC_Binary_HW_Share<float>;

}

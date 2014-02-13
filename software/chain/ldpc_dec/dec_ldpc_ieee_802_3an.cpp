//
//  Copyright (C) 2012 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Hardware-compliant LDPC decoder module for the IEEE 802.3an standard
/// \author Matthias Alles
/// \date   2012/03/26
//

#include "dec_ldpc_ieee_802_3an.h"
#include "dec_ldpc_ieee_802_3an_codes.hpp"

using namespace hlp_fct::logging;
using namespace cse_lib::ieee_802_3an_codes;
using namespace std;

namespace cse_lib {

template<class T>
void Decoder_LDPC_IEEE_802_3an<T>::Set_LDPC_Parameters()
{

	/*
	 * Parameterize the _Share class with the parameters from the configuration.
	 */

	check_node_algorithm_ = dec_algorithm();
	num_lambda_min_       = num_lambda_min();
	esf_factor_           = esf_factor();
	bw_fract_             = bw_fract();
	num_partitions_       = num_partitions();
	threshold_            = threshold();

	// Calculate the maximum values that can be represented by the chosen quantization.
	max_msg_extr_         = (1 << (bw_extr() - 1)) - 1;
	max_msg_app_          = (1 << (bw_app()  - 1)) - 1;

	num_variable_nodes_ = 2048;
	num_check_nodes_ = 384;
	src_parallelism_ = 1;
	dst_parallelism_ = 1;
	max_check_degree_ = 32;
	is_IRA_code_ = false;
	addr_vector_ = &ieee_802_3an_p1_n2048_r084_addr[0];
	shft_vector_ = &ieee_802_3an_p1_n2048_r084_shft[0];
}


template<class T>
void Decoder_LDPC_IEEE_802_3an<T>::Init()
{

	// Set code and decoder parameters
	Set_LDPC_Parameters();

	//mean_iterations.Reset_Port();
	mean_iterations.Reset();
	flipped_bits.Reset();

	// Resize output buffers and internal RAMs.
	try
	{
		// Output RAM content for each iteration.
		output_bits().Resize(num_iterations(), num_variable_nodes_);
		output_bits_llr_app().Resize(num_iterations(), num_variable_nodes_);

		// Decoder RAMs
		app_ram_.Resize(dst_parallelism_, num_variable_nodes_ / dst_parallelism_);
		msg_ram_.Resize(dst_parallelism_, num_check_nodes_ * max_check_degree_ / dst_parallelism_);
	}
	catch(runtime_error& e)
	{
//		string str = Exception_Str_Memory_Alloc("Decoder_LDPC_IEEE_802_3an<T>", instance_name(), e.what());
//		throw runtime_error(str);
	}

	param_list_.config_modified(false);
	input_data_list_.port_modified(false);
}


template<class T>
int Decoder_LDPC_IEEE_802_3an<T>::Run()
{

	unsigned int pchk_satisfied;
	unsigned int iter = 0;
	bool next_iter_is_last_iter = false;
	bool last_iter = false;

	decoding_successful().Write(false);
	num_modified_systematic_bits().Write(0);

	if(param_list_.config_modified())
		Init();

	// Read the channel values and store them in app_ram_.
	this->Init_APP_RAM(is_IRA_code_, input_bits_llr(), app_ram_);

	do
	{
		// Perform one ldpc decoder iteration.
		switch(scheduling())
		{
		case LAYERED:
			pchk_satisfied = this->Decode_Layered(app_ram_, msg_ram_, iter);
			break;

		case TWO_PHASE:
			pchk_satisfied = this->Decode_Two_Phase(app_ram_, msg_ram_, iter, app_parity_check());
			break;

		default:
			pchk_satisfied = 0;
			Msg(ERROR, instance_name(), "Selected scheduling not supported for these codes!");
			break;
		}

		/*
		 * Read the app_ram_ and store APP values in output_bits_llr_app() and
		 * hard decoded bits in output_bits buffer.
		 */
		this->Read_APP_RAM(app_ram_, iter, output_bits_llr_app(), output_bits());

		// Check whether all parity checks were satisfied in the previous iteration.
		last_iter = next_iter_is_last_iter;

		// Are all parity checks satisfied?
		if (pchk_satisfied == num_check_nodes_)
		{
			decoding_successful().Write(true);
			next_iter_is_last_iter = true; // Do one more iteration, as it is done in hardware!
		}

 	// Store the number of flipped bits
        if (iter != 0) {
            flipped_bits(iter)().Write(this->Calc_Flipped_Bits(iter, output_bits()));
        }

		// Increase iteration counter.
		iter++;

		mean_iterations(iter)().Write(iter);


	/*
	 * Abort conditions:
	 * 1) maximum number of iterations is reached
	 * 2) all parity checks are satisfied: Since the hardware performs one more
	 * iteration after all parity checks are satisfied, we delay the stopping
	 * in the software as well.
	 */
	} while (iter < num_iterations() &&
	         last_iter == false);

	// Write the number of unsatisfied parity checks.
	num_unsatisfied_parity_checks().Write(num_check_nodes_ - pchk_satisfied);

	// Set number of used iterations in output buffer.
	iterations_performed().Write(iter);

	// Get statistic about modified bits.
	num_modified_systematic_bits().Write(this->Calc_Modified_Systematic_Bits(iter, input_bits_llr(), output_bits()));

	// Fill the output buffer and the status port for the remaining iterations.
	for(unsigned int i = iter; i < num_iterations(); i++)
	{
		mean_iterations(i + 1)().Write(iter);
		output_bits_llr_app()[i] = output_bits_llr_app()[iter - 1];
		output_bits()[i]         = output_bits()[iter - 1];
	}

	return 0;
}

template class Decoder_LDPC_IEEE_802_3an<int>;
template class Decoder_LDPC_IEEE_802_3an<float>;
}

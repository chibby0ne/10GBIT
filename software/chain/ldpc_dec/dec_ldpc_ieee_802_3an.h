//
//  Copyright (C) 2010 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Hardware-compliant LDPC decoder module for the IEEE 802.15.3c standard.
/// \author Matthias Alles
/// \date   2010/12/13
//

#ifndef DEC_LDPC_IEEE_802_3AN_H_
#define DEC_LDPC_IEEE_802_3AN_H_

#include "dec_ldpc_bin_hw_share.h"
#include "dec_ldpc_bin_hw_iface.h"
#include "dec_ldpc_ieee_802_3an_param.h"


namespace cse_lib {

/// LDPC decoder module for the IEEE_802_3AN Standard (N = 2048, P = 1, R = 0.84 ).
/**
 * \ingroup modules
 */
template<typename T>
class Decoder_LDPC_IEEE_802_3an : public Decoder_LDPC_Binary_HW_Interface<T>,
                                  public Decoder_LDPC_IEEE_802_3an_Parameter,
                                  public Decoder_LDPC_Binary_HW_Share<T>
{

public:

	Decoder_LDPC_IEEE_802_3an() { };
	virtual ~Decoder_LDPC_IEEE_802_3an() { };

	int Run();

	// Make members of template classes visible.
	using Decoder_LDPC_Binary_HW_Interface<T>::input_data_list_;
	using Decoder_LDPC_Binary_HW_Interface<T>::input_bits_llr;
	using Decoder_LDPC_Binary_HW_Interface<T>::output_bits;
	using Decoder_LDPC_Binary_HW_Interface<T>::output_bits_llr_app;
	using Decoder_LDPC_Binary_HW_Interface<T>::iterations_performed;
	using Decoder_LDPC_Binary_HW_Interface<T>::flipped_bits;
	using Decoder_LDPC_Binary_HW_Interface<T>::mean_iterations;
	using Decoder_LDPC_Binary_HW_Interface<T>::decoding_successful;
	using Decoder_LDPC_Binary_HW_Interface<T>::num_modified_systematic_bits;
	using Decoder_LDPC_Binary_HW_Interface<T>::num_unsatisfied_parity_checks;

private:

	Buffer<T, 2> app_ram_;  ///< APP RAM of LDPC decoder
	Buffer<T, 2> msg_ram_;  ///< Extrinsic RAM of LDPC Decoder

	void Init();


	/// Set decoder and code parameters
	/**
	 * Sets default decoder parameters to:
	 * - Layered decoding
	 * - Min-Sum algorithm with extrinsic scaling factor (0.875)
	 * - 5 iterations
	 * - 6 bits for extrinsics
	 * - 6 bits for channel values
	 * - 8 bits for APP values
	 */
	void Set_LDPC_Parameters();

	// Make members of template classes visible.
	using Decoder_LDPC_Binary_HW_Share<T>::check_node_algorithm_;
	using Decoder_LDPC_Binary_HW_Share<T>::num_lambda_min_;
	using Decoder_LDPC_Binary_HW_Share<T>::esf_factor_;
	using Decoder_LDPC_Binary_HW_Share<T>::num_partitions_;
	using Decoder_LDPC_Binary_HW_Share<T>::threshold_;
	using Decoder_LDPC_Binary_HW_Share<T>::bw_fract_;
	using Decoder_LDPC_Binary_HW_Share<T>::max_msg_extr_;
	using Decoder_LDPC_Binary_HW_Share<T>::max_msg_app_;
	using Decoder_LDPC_Binary_HW_Share<T>::num_variable_nodes_;
	using Decoder_LDPC_Binary_HW_Share<T>::num_check_nodes_;
	using Decoder_LDPC_Binary_HW_Share<T>::max_check_degree_;
	using Decoder_LDPC_Binary_HW_Share<T>::src_parallelism_;
	using Decoder_LDPC_Binary_HW_Share<T>::dst_parallelism_;
	using Decoder_LDPC_Binary_HW_Share<T>::is_IRA_code_;
	using Decoder_LDPC_Binary_HW_Share<T>::addr_vector_;
	using Decoder_LDPC_Binary_HW_Share<T>::shft_vector_;
};
}
#endif // DEC_LDPC_IEEE_802_3AN_H_

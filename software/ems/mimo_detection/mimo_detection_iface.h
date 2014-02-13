//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * mimo_detection_iface.h
 *
 *  Created on: 03.02.2011
 *      Author: gimmler
 */

#ifndef MIMO_DETECTION_IFACE_H_
#define MIMO_DETECTION_IFACE_H_

#include "../modules/base/base_iface.h"

using namespace cse_lib;

///interface for all mimo detection algorithms
template<class T>
class Mimo_Detection_Interface: public Base_Interface {
public:
	Mimo_Detection_Interface() {
		//register inputs and outputs
		input_R_matrices.Register("input_R_matrices", input_data_list_);
		input_symb_yhat.Register("input_symbols_y_hat", input_data_list_);
		input_permutation_vector.Register("input_permutation_vector", input_data_list_);
		input_bits_llr_apriori.Register("input_bits_llr_apriori", input_data_list_);

		output_bits_llr_app.Register("output_bits_llr_app", output_data_list_);
		output_bits_llr_extrinsic.Register("output_bits_llr_extrinsic", output_data_list_);
		output_bits.Register("output_bits", output_data_list_);
		output_bits_extr.Register("output_bits_extr", output_data_list_);

		mean_num_bm_calcs.Register("mean_num_bm_calcs", status_out_list_);
		mean_num_leaves.Register("mean_num_leaves", status_out_list_);
		mean_num_nodes.Register("mean_num_nodes", status_out_list_);
		hist_num_bm_calcs.Register("hist_num_bm_calcs", status_out_list_, false);
		hist_num_leaves.Register("hist_num_leaves", status_out_list_, false);
		hist_num_nodes.Register("hist_num_nodes", status_out_list_, false);

	}
	virtual ~Mimo_Detection_Interface() {
	}

	///Input R matrices
	Data_In<complex<T> , 2> input_R_matrices;
	///Input transformed channel values
	Data_In<complex<T> > input_symb_yhat;
	///Input permutation vectors for each transmitted vector PermutationVector(index_y) = index_y_p
	Data_In<int, 2> input_permutation_vector;
	///input apriori llrs
	Data_In<T> input_bits_llr_apriori;

	///output app llrs
	Data_Out<T> output_bits_llr_app;
	///output extrinsic llrs
	Data_Out<T> output_bits_llr_extrinsic;
	///output hard bits (for error rates only)
	Data_Out<unsigned int, 2> output_bits;
	///output hard bits (for error rates only)
	Data_Out<unsigned int, 2> output_bits_extr;

	//Stat_Out ports
	///average number of visited nodes
	Status_Out<long long int, 1, Status_Out_Plugin_Mean> mean_num_nodes;
	///average number of branch metric calculations
	Status_Out<long long int, 1, Status_Out_Plugin_Mean> mean_num_bm_calcs;
	///average number of leaf nodes
	Status_Out<long long int, 1, Status_Out_Plugin_Mean> mean_num_leaves;

	///histogram of visited nodes
	Status_Out<int, 1, Status_Out_Plugin_Histogram> hist_num_nodes;
	///histogram of branch metric calculations
	Status_Out<int, 1, Status_Out_Plugin_Histogram> hist_num_bm_calcs;
	///histogram of leaf nodes
	Status_Out<int, 1, Status_Out_Plugin_Histogram> hist_num_leaves;

};

#endif /* MIMO_DETECTION_IFACE_H_ */

//
//  Copyright (C) 2010 - 2013 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Hardware-compliant LDPC decoder shared functionality.
/// \author Matthias Alles
/// \date   2010/12/03
//
//
#ifndef DEC_LDPC_BIN_HW_SHARE_H_
#define DEC_LDPC_BIN_HW_SHARE_H_

#include <cstdlib>
#include "../../cse/assistance/buffer.h"
#include "dec_ldpc_check_node_share.h"
#include <vector>


namespace cse_lib {

/// This is shared functionality for hardware compliant LDPC decoding.
/**
 * \ingroup share
 */
template<typename T>
class Decoder_LDPC_Binary_HW_Share : public Decoder_LDPC_Check_Node_Share
{

public:

	Decoder_LDPC_Binary_HW_Share() { };
	virtual ~Decoder_LDPC_Binary_HW_Share() { };

protected:


	/**************************************************************************************
	 ** Code describing values, have to be set before using the functions of this class! **
	 **************************************************************************************/

	/*
	 * Decoding parameters
	 */

	/// Check node algorithm
	CHECK_NODE_ENUM check_node_algorithm_;

	/// Number of messages to take in account for lambda-min algorithm
	unsigned int num_lambda_min_;

	/// Extrinsic scaling factor for Min-Sum algorithm.
	/**
	 * For fix-point implementations we support 0.875 and 0.75.
	 * For float implementations any value is valid.
	 */
	float esf_factor_;

  /// Number of Partitions (Split Row)
    unsigned int num_partitions_;

    /// Threshold (threshold split Row)
    unsigned int threshold_;



	/*
	 * Quantization
	 */

	/// Number of bits for fractional part of all values. Relevant for some fix-point check node algorithms.
	unsigned int bw_fract_;

	/// Maximum value of extrinsic information, only relevant for fix-point implementations.
	unsigned int max_msg_extr_;

	/// Maximum value of APP values, only relevant for fix-point implementations.
	unsigned int max_msg_app_;


	/*
	 * Code Parameters
	 */

	/// Number of variable nodes of the code (= N)
	unsigned int num_variable_nodes_;

	/// Number of check nodes of the code (= N - K = M)
	unsigned int num_check_nodes_;

	/// Maximum check node degree of the code
	unsigned int max_check_degree_;

	/// Parallelism of the base matrix as defined in the standard, i.e., submatrix size
	unsigned int src_parallelism_;

	/// Currently implemented parallelism of the code, i.e., actual decoder parallelism
	unsigned int dst_parallelism_;

	/// IRA codes have only parity nodes of degree two and one parity node of degree one.
	bool is_IRA_code_;

	/// Pointer to the address vector of the currently processed code
	int * addr_vector_;

	/// Pointer to the shift vector of the currently processed code
	int * shft_vector_;


	/****************************************************
	 ** Functions for hardware compliant LDPC decoding **
	 ****************************************************/


	/// Modify left shifts to right shifts, in order to be compliant with the C++ code.
	/**
	 * The decoder functions of this class only perform right shifting of
	 * submatrices.  Some standards, however, define the submatrices with left
	 * shifts instead. By calling this function, before processing the left
	 * shifts as defined in the standard are modified to be right shifts, as
	 * required by the decoder functions.
	 */
	void Change_Shift_Direction();


	/// Check node implementation according to Lambda-Min algorithm
	/**
	 * This function takes quantized values and performs the check node
	 * operation according to the Lambda-Min approximation.
	 *
	 * \param in_out_msg     Input/Output array of messages to/from the check
	 *                       node.
	 * \param num_lambda_min Lambda, either 2 or 3.
	 * \param num_bits_fract Number of fractional bits of input values.
	 *
	 * \return -1 : Parity check was not satisfied.
	 * \return  1 : Parity check was satisfied.
	 */
	int Check_Node_Lambda_Min(Buffer<T> &in_out_msg,
	                          unsigned int num_lambda_min,
	                          unsigned int num_bits_fract);


	/// Check node implementation according to Min-Sum algorithm
	/**
	 * This function takes quantized values and performs the check node
	 * operation according to the Min-Sum approximation. The extrinsic scaling
	 * factor esf_factor is used for output scaling.
	 *
	 * \param in_out_msg   Input/Output array of messages to/from the check node.
	 * \param esf_factor   Extrinsic scaling factor: Supported are 0.75 and 0.875.
	 *
	 * \return -1 : Parity check was not satisfied.
	 * \return  1 : Parity check was satisfied.
	 */
	int Check_Node_Min_Sum(Buffer<T> &in_out_msg, float esf_factor);

 /// Check node implementation according to Multi-Split-Row
	/**
	 * This function takes quantized values and performs the check node
	 * operation according to the Min-Sum approximation. The extrinsic scaling
	 * factor esf_factor is used for output scaling.
	 *
	 * \param in_out_msg[]   Input/Output array of messages to/from the check
	 *                       node.
	 * \param esf_factor     Extrinsic scaling factor: Supported are 0.75 and 0.875.
     *
     * \param num_partitions Number of partitions to divide Parity Check Matrix
	 *
	 * \return -1 : Parity check was not satisfied.
	 * \return  1 : Parity check was satisfied.
	 */
	int Check_Node_Split_Row(Buffer<T> &in_out_msg, float esf_factor, unsigned int num_partitions, int threshold);
    
    void Distribute_Messages_Partitions(std::vector< std::vector<T> > &parts, Buffer<T> &in_out_msg, unsigned int partitions, unsigned int len);

    int Split_Row_Local_Minimum_Sign(std::vector<T>& part, T& minimum1, unsigned int& min1_idx, T& minimum2, unsigned int& threshold_en, int threshold, unsigned int part_num);

    bool Neighbors_Threshold_Enable(unsigned int part_num, unsigned int partitions, std::vector<unsigned int>& treshold_en);
//    bool Neighbors_Threshold_Enable(unsigned int part_num, unsigned int partitions, unsigned int treshold_en[]);

    void Update_Minimum(std::vector<T>& temp_part, unsigned int min1_idx, T min1, T min2, unsigned int len);

    void Multiply_Sign_Esf(std::vector<T>& part, std::vector<T>& temp_part, int local_pchk_sign, float esf_factor, unsigned int part_num);

    void Multiply_Split_Sign_Signal(std::vector<T>& part, int split_sign);

    void Regroup_Messages_Partitions(std::vector< std::vector<T> > &parts, Buffer<T> &in_out_msg, unsigned int len);

    unsigned int Minimum_Fixp_Split_Row(Buffer<T>& message, unsigned int& idx, unsigned int part_num, unsigned int len);

    unsigned int Calc_Flipped_Bits(unsigned int iter, Buffer<unsigned int, 2> &output_bits);





	/// Check node implementation according to the 3-Min algorithm with optimized data path.
	/**
	 * This function takes quantized values and performs the check node
	 * operation according to the Lambda-3Min approximation. The output values are generated from
	 * 2 of the 3 values only. This hardware implementation friendly algorithm loses a little bit
	 * communications performance compared to the original Lambda 3Min implementation.
	 *
	 * \param in_out_msg       Input/Output array of messages to/from the check node.
	 * \param num_bits_fract   Number of fractional bits of input values.
	 *
	 * \return -1 : Parity check was not satisfied.
	 * \return  1 : Parity check was satisfied.
	 */
	int Check_Node_Lambda_3Min_3Mag(Buffer<T> &in_out_msg, unsigned int num_bits_fract);

	/// Copy input data to LDPC decoder RAMs.
	/**
	 * Function that initializes the int app_ram_ with the quantized values stored
	 * in input_bits_llr(). The operation performs the interleaving that is
	 * dependant on the decoder parallelism. This function has to be called once
	 * before the decoder function decode_xxx() can be called for the first time.
	 *
	 * \param parity_reordering  Perform parity reordering as required by IRA codes.
	 * \param input_bits_llr     LLR values of the codeword in linear order.
	 * \param app_ram            Buffer to store LLR values as in hardware decoder.
	 */
	void Init_APP_RAM(bool         parity_reordering,
	                  Buffer<T>   &input_bits_llr,
	                  Buffer<T, 2> &app_ram);


	/// Copy LDPC decoder RAMs deinterleaved to output arrays.
	/**
	 * This function has to be called when the soft information or the hard decoded
	 * bits after decoding have to be read. This can happen after each decoder
	 * iteration. If only the result after the maximum number of iterations is of
	 * interested calling only once is sufficient.
	 * Deinterleaving according to the decoder parallelism is performed and either
	 * a 0 or a 1 is stored for the hard decoded bits.
	 * Modifies output_bits() and output_bits_llr_app().
	 *
	 * \param[in]   app_ram  Two-dimensional Buffer containing APP values after decoding.
	 * \param[in]   iter     Current LDPC decoder iteration
	 * \param[out] output_bits_llr_app data container to store the soft output information.
	 * \param[out] output_bits data container to store the hard output information.
	 */
	void Read_APP_RAM(Buffer<T, 2>   &app_ram,
	                  int             iter,
	                  Buffer<T, 2>   &output_bits_llr_app,
	                  Buffer<unsigned int, 2> &output_bits);


	/// Calculate the number of systematic bits that have been swapped during decoding.
	/**
	 *
	 * Compares input and output buffers of the decoder and calculates how many bits
	 * have been swapped during decoding (hard).
	 *
	 * \param iter            The iteration to check for.
	 * \param input_bits_llr  The input llr values of the LDPC decoder, linear order.
	 * \param output_bits     The hard decoded bits of the LDPC decoder, for each iteration.
	 */
	unsigned int Calc_Modified_Systematic_Bits(unsigned int             iter,
	                                           Buffer<T>               &input_bits_llr,
	                                           Buffer<unsigned int, 2> &output_bits);


	/// Layered LDPC decoder function. Each call corresponds to a single iteration.
	/**
	 * This function contains the LDPC decoder itself. Each call corresponds to a
	 * single decoder iteration. So calling this function in a loop will result
	 * in the normal LDPC decoding process. The memories can be accessed after
	 * each iteration. The decoder returns the number of parity checks that were
	 * satisfied which allows to stop the decoding process when all parity checks
	 * are satisfied.
	 *
	 * This particular function corresponds to hardware-compliant decoding in
	 * case of superposed submatrices in the parity check matrix. These exist
	 * for instance in DVB-S2/T2/C2 and the CMMB standard. Superposed submatrix
	 * means, a single check node group accesses more than once to the same
	 * variable node group, i.e., the address vector for one check node group
	 * contains the same entry more than once.
	 *
	 * \param app_ram  Buffer of two dimensional memory array storing the
	 *                 quantized APP LLR values.
	 * \param msg_ram  Buffer of two dimensional memory array storing the
	 *                 extrinsic information (can be uninitialized for the
	 *                 first iteration).
	 * \param iter     The iteration the decoder currently processes.
	 *
	 * \return Number of parity checks that were satisfied.
	 */
	unsigned int Decode_Layered_Superposed(Buffer<T, 2> &app_ram,
	                                       Buffer<T, 2> &msg_ram,
	                                       int             iter);


	/// Layered LDPC decoder function. Each call corresponds to a single iteration.
	/**
	 * This function contains the LDPC decoder itself. Each call corresponds to a
	 * single decoder iteration. So calling this function in a loop will result
	 * in the normal LDPC decoding process. The memories can be accessed after
	 * each iteration. The decoder returns the number of parity checks that were
	 * satisfied which allows to stop the decoding process when all parity checks
	 * are satisfied.
	 *
	 * This particular function corresponds to hardware-compliant decoding with
	 * a layered decoding schedule.
	 *
	 * \param app_ram  Buffer of two dimensional memory array storing the
	 *                 quantized APP LLR values.
	 * \param msg_ram  Buffer of two dimensional memory array storing the
	 *                 extrinsic information (can be uninitialized for the
	 *                 first iteration).
	 * \param iter     The iteration the decoder currently processes.
	 *
	 * \return Number of parity checks that were satisfied.
	 */
	unsigned int Decode_Layered(Buffer<T, 2> &app_ram,
	                            Buffer<T, 2> &msg_ram,
	                            int             iter);


	/// Two-Phase LDPC decoder function. Each call corresponds to a single iteration.
	/**
	 * This function contains the LDPC decoder itself. Each call corresponds to a
	 * single decoder iteration. So calling this function in a loop will result
	 * in the normal LDPC decoding process. The memories can be accessed after
	 * each iteration. The decoder returns the number of parity checks that were
	 * satisfied which allows to stop the decoding process when all parity checks
	 * are satisfied.
	 *
	 * This particular function corresponds to hardware-compliant decoding with
	 * a two-phase scheduling. All check nodes work on the same data and no
	 * early update occurs.
	 *
	 * \param app_ram  Buffer of two dimensional memory array storing the
	 *                 quantized APP LLR values.
	 * \param msg_ram  Buffer of two dimensional memory array storing the
	 *                 extrinsic information (can be uninitialized for the
	 *                 first iteration).
	 * \param iter     The iteration the decoder currently processes.
	 * \param app_parity_check Defines whether to use APP values or extrinsic
	 *                         values for calculation of the parity check.
	 *
	 * \return Number of parity checks that were satisfied.
	 */
	unsigned int Decode_Two_Phase(Buffer<T, 2> &app_ram,
	                              Buffer<T, 2> &msg_ram,
	                              int             iter,
	                              bool            app_parity_check = false);


	/// Function retrieve extrinsic values to feed the chcek nodes with.
	/**
	 * \param app_ram        Buffer to read APP values from
	 * \param msg_ram        Extrinsic storage to read check node data from previous iteration
	 * \param iter           Current decoder iteration (starting from 0)
	 * \param cng_counter    Check node group counter
	 * \param cfu_counter    Check node functional unit counter (within the check node group)
	 * \param check_node_in  Values that are going to the check node
	 * \param app_out        APP values that were used to calculate the check node input data
	 */
	void Get_Check_Node_Input(Buffer<T, 2> &app_ram,
	                          Buffer<T, 2> &msg_ram,
	                          unsigned int iter,
	                          unsigned int cng_counter,
	                          unsigned int cfu_counter,
	                          Buffer<T>  &check_node_in,
	                          Buffer<T>  &app_out);


	/// Function to write back APP values to app_ram.
	/**
	 * \param app_ram        Buffer to write APP values to
	 * \param msg_ram        Extrinsic storage to write check node data into
	 * \param iter           Current decoder iteration (starting from 0)
	 * \param cng_counter    Check node group counter
	 * \param cfu_counter    Check node functional unit counter (within the check node group)
	 * \param check_node_out Buffer holding the information of the check node we processed
	 */
	void Write_Check_Node_Output(Buffer<T, 2> &app_ram,
	                             Buffer<T, 2> &msg_ram,
	                             unsigned int iter,
	                             unsigned int cng_counter,
	                             unsigned int cfu_counter,
	                             Buffer<T>  &check_node_out);


	/// Calculate parity check in input buffer
	/**
	 * \param values  Buffer to calculate parity check on.
	 *
	 * \return  1 : Parity check was not satisfied.
	 * \return  0 : Parity check was satisfied.
	 */
	inline unsigned int Calc_Parity_Check(Buffer<T> &values);


private:


	/// Search for minimum in an array.
	/**
	 * \param message           Buffer of input values.
	 * \param idx               Position of the minimum in the array.
	 *
	 * \return Returns the minimum within the array and the position (idx).
	 */
	inline T Minimum_Fixp(Buffer<T> &message, unsigned int &idx);


	/// Calculate the delta function in a hardware compliant version.
	/**
	 * Calulates
	 * \f[
	 * \ln \left( 1 + e^{-|x|}\right) \approx max \left( \frac{5}{8} - \frac{|x|}{4}, 0 \right)
	 * \f]
	 *
	 * \param x               Input value, has to be positive.
	 * \param num_bits_fract  Number of fractional bits of input values.
	 *
	 * \return Return value is computed according to delta function (see above).
	 */
	inline int Delta(int x, unsigned int num_bits_fract);


	/// Calculates the delta function in a float version.
	/**
	 * Calulates
	 * \f[
	 * \ln \left( 1 + e^{-|x|}\right)
	 * \f]
	 *
	 * \param x               Input value, has to be positive.
	 * \param num_bits_fract  Unused.
	 *
	 * \return Return value is computed according to delta function (see above).
	 */
	inline float Delta(float x, unsigned int);


	/// Calculate Min* of two input values.
	/**
	 * This function looks for the minimum of a and b and adds a correction
	 * term to the minimum. A fixpoint representation is assumed where frac_bits
	 * gives the number of fractional bits within the integer input values
	 *
	 * \f[
	 *     min^*(a,b) = min(a,b) + \ln \left( 1 + e^{-|a+b|}\right)
	 *                           - \ln \left( 1 + e^{-|a-b|}\right)
	 * \f]
	 *
	 * In hardware the logarithmus-terms are approximated by the \f$\delta\f$
	 * function, see * delta().
	 *
	 * \param a              Input value 1.
	 * \param b              Input value 2.
	 * \param num_bits_fract Number of fractional bits of input values.
	 *
	 * \return Minumum of a and b plus correction term.
	 */
	inline int Minstar_Fixp(T a, T b, unsigned int num_bits_fract);

	Buffer<T>   cn_msg_abs_;
	Buffer<int> cn_msg_sign_;

	/// Scales with 0.875 or 0.75 as it is done in hardware (bit-accurate!)
	inline void Scale_ESF(int &msg, float esf_factor);

	/// Simple Multiplication of msg with esf_factor.
	inline void Scale_ESF(float &msg, float esf_factor);

	/// Saturate symmetrically to be hardware compliant.
	inline void Saturate(int &input, int max_msg);

	/// For float models, we do not saturate.
	inline void Saturate(float& , float);

	/// Saturate symmetrically to be hardware compliant.
	inline void Saturate(int &input, int min_msg, int max_msg);

	/// For float models, we do not saturate.
	inline void Saturate(float&, float, float);

};
}
#endif // DEC_LDPC_BIN_HW_SHARE_H_

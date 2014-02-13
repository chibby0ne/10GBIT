//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * conv_share.h
 *
 *  Created on: 15.02.2011
 *      Author: gimmler
 */

#ifndef CONV_SHARE_H_
#define CONV_SHARE_H_

#include <math.h>
#include "itpp/base/mat.h"
#include "itpp/base/vec.h"
#include "itpp/base/converters.h"
#include "../assistance/buffer.h"

using namespace itpp;
using namespace std;

///shared functionalities for convolutional codes
/**
 * This class computes transition tables depending on two forward polynomials and one feedback polynomial.
 * It can also compute the number of states or M.
 *
 * ATTENTION: currently this class is not generic!!! It only returns values for a 64-state NSC with polynomials 133 and 171.
 */
class Convolutional_Share {
public:

protected:
	Convolutional_Share() {
	}
	virtual ~Convolutional_Share() {
	}

	///state_transition_table[current_state][input_bit] contains new_state
	Buffer<unsigned int, 2> state_transition_table;
	///out_par1_table[current_state][input_bit] contains the resulting value for parity1
	Buffer<unsigned int, 2> out_par1_table;
	///out_par1_table[current_state][input_bit] contains the resulting value for parity2
	Buffer<unsigned int, 2> out_par2_table;

	///computes the tables state_transition_table, out_par1_table, out_par2_table2 depending on the given polynomials
	void Get_tables(unsigned int forward_polynomial1, unsigned int forward_polynomial2, unsigned int feedback_polynomial = 0);
	///computes the number of states depending on the given polynomials
	unsigned int Get_num_states(unsigned int forward_polynomial1, unsigned int forward_polynomial2, unsigned int feedback_polynomial);
	///computes M depending on the given polynomials
	int Get_M(unsigned int forward_polynomial1, unsigned int forward_polynomial2, unsigned int feedback_polynomial);

private:
	///transforms a decimal number i into a binary number with stellen bits, returns the binary number as bvec
	bvec transform_dec_to_dual(int i, int stellen);
	///transform dual number in into an integer
	int transform_dual_to_dec(bvec in);
};

#endif /* CONV_SHARE_H_ */

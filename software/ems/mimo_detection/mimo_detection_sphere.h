//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * mimo_detection_sphere.h
 *
 *  Created on: 03.02.2011
 *      Author: gimmler
 */

#ifndef MIMO_DETECTION_SPHERE_H_
#define MIMO_DETECTION_SPHERE_H_

#include "itpp/itbase.h"
#include "itpp/itcomm.h"

#include "../assistance/buffer.h"

#include "mimo_detection_iface.h"
#include "mimo_detection_sphere_param.h"

using namespace itpp;
using namespace cse_lib;

///flexible mimo sphere decoder
template<class T>
class Mimo_Detection_Sphere: public Mimo_Detection_Interface<T> , public Mimo_Detection_Sphere_Parameter, public Mapper_Share {
public:
	Mimo_Detection_Sphere();
	virtual ~Mimo_Detection_Sphere();

	///process one block of input data
	int Run();
	///process one block of input data
	int Run(int Iteration);

	///iteration in which channel decoder starts after last outer iteration
	int idx_final_inner_iterations;

protected:
	///compute derived parameters and resize memories
	void Init();
	///resize all memories
	void Alloc_Mem();

//private:


	///number of modulation symbols per block
	int SymbolsPerBlock;
	///number of bits per modulation symbol
	unsigned int BitsPerSymbol;
	///number of channel uses per block
	int NrTransmissions;

	///counts the number of leaves during one detection
	unsigned int leaves_ctr;
	///counts the number of visited nodes during on detection
	unsigned int nodes_ctr;
	///counts the number of computed branch metrics during on detection
	unsigned int bm_calcs_ctr;
	///tree level dependent sphere constraint for APP_LEVEL
	vec AppPerLevel;
	///list of best |K|=ListSize solutions found so far
	Buffer<T> ListOfBestSolutions;
	///apriori information for each possible modulation symbol in each level of the tree
	Buffer<T, 2> tree_apriori_per_symbol;
	///maximum apriori symbol information for each layer
	Buffer<T> max_apriori_symbol_per_layer;
	///enumeration technique per layer
	ENUMERATION_TECHNIQUE *enum_technique_per_layer;
	///enumeration order for apriori enumeration
	Buffer<T,2> apriori_enumeration;
	///apriori LLRs for tree search
	T *tree_apriori_llrs;
	///APP LLRs, result of tree search
	T *tree_app_llrs;
	///channel values
	complex<T> *YHatTree;

	///Table of all modulation symbols as bit vectors
	Buffer<bool, 2> ConstellationBits;
	///Modulation symbols according to ConstellationBits
	Buffer<complex<T> > ConstellationPoints;
	///constants for bias reduction
	Buffer<T> sigma_square_si_square_table_;
	///precalculated multiplications of matrix elements and modulation symbols rij_sk_table[i][j][k] = r_ij * s_k
	complex<T> ***rij_sk_table_;
	///allocate a 3 dimensional array with dimensions dim1, dim2, dim3
	template<class U> U *** allocate3Darray(unsigned int dim1, unsigned int dim2, unsigned int dim3) {
		U ***p3DArray;
		// Allocate memory
		p3DArray = new U**[dim1];
		for (unsigned int i = 0; i < dim1; ++i) {
			p3DArray[i] = new U*[dim2];
			for (unsigned int j = 0; j < dim2; ++j) {
				p3DArray[i][j] = new U[dim3];
			}
		}
		//set all entries to zero
		for (unsigned int i = 0; i < dim1; i++)
			for (unsigned int j = 0; j < dim2; j++)
				for (unsigned int k = 0; k < dim3; k++)
					p3DArray[i][j][k] = 0;
		return p3DArray;
	}
	///allocate a 2 dimensional array with dimensions dim1, dim2
	template<class U> U ** allocate2Darray(unsigned int dim1, unsigned int dim2) {
		U **p2DArray;
		// Allocate memory
		p2DArray = new U*[dim1];
		for (unsigned int i = 0; i < dim1; ++i) {
			p2DArray[i] = new U[dim2];
		}
		//set all entries to zero
		for (unsigned int i = 0; i < dim1; i++)
			for (unsigned int j = 0; j < dim2; j++)
				p2DArray[i][j] = 0;
		return p2DArray;
	}
	///delete a 3 dimensional array with dimensions dim1, dim2, dim3
	template<class U> void deallocate3Darray(U ***p3DArray, unsigned int dim1, unsigned int dim2, unsigned int dim3) {
		for (unsigned int i = 0; i < dim1; ++i) {
			for (unsigned int j = 0; j < dim2; ++j) {
				delete[] p3DArray[i][j];
			}

			delete[] p3DArray[i];
		}
		delete[] p3DArray;
	}
	///delete a 2 dimensional array with dimensions dim1, dim2
	template<class U> void deallocate2Darray(U **p2DArray, unsigned int dim1, unsigned int dim2) {
		for (unsigned int i = 0; i < dim1; ++i) {
			delete[] p2DArray[i];
		}
		delete[] p2DArray;
	}

	///precalculates all possible values rij * sj for one matrix R
	inline void precalculate_rij_sk_table(complex<T> ***rij_sk_tab, complex<T> **matrix_r);
	///precalculates the table sigma_square * symbol_energy
	inline void precalculate_sigma_square_si_square_table(Buffer<T> &sigma_square_si_square_table, T sigmaa_square);
	///precalculates the apriori information per symbol
	inline void precalculate_apriori_per_symbol(Buffer<T, 2> &apriori_per_symbol, T *apriori_per_bit, Buffer<T> &sigma_square_si_square_table);

	///generates the constellation tables ConstellationBits and ConstellationPoints, called by constructor
	virtual inline void generate_constellation(Mapper_Share::MAPPINGS_ENUM, Buffer<bool, 2> &constellation_bits, Buffer<complex<T> > &constellation_points);
	/**
	 * \brief computes LLRPermuteIndex out of inPermutationVector
	 *
	 * @param t index of transmission vector
	 */
	inline void CalcLLRPermuteIndex(int *permutation_vector, int *permutation_bit_index);
	///clip llrs to output_clipping_value
	inline void ClipLlrs(T *LlrVec, int length);
	/**
	 * \brief computes apriori information for one symbol
	 *
	 * @param apriori_llrs apriori LLRs for the according bits
	 * @param test_bits bit vector representing modulation symbol
	 * @return apriori information of symbol
	 */
	inline T apriori_update(T *apriori_llrs, bool *test_bits);
	/**
	 * \brief fully enumerated tree search
	 *
	 * @param t index of transmission vector
	 */
	void tree_search_depth_first(int t);

	void tree_search_sic(int t);
	/**
	 * \brief depth-first tree search without enumeration
	 *
	 * @param t index of transmission vector
	 */
	void tree_search_without_enumeration(int t);
	/**
	 * \brief calculation of partial euclidian distance (PED)
	 *
	 * @param p path in tree
	 * @param d level in tree
	 * @param PED partial euclidian distance
	 * @param b_i interference reduced channel value
	 * @param row_offset offset value to read from inPreCalculatedr_s
	 */
	inline void get_PED(unsigned int *p, int d, T *PED, complex<T> b_i);
	/**
	 * \brief calculation of interference reduced channel value
	 *
	 * @param p path in tree
	 * @param d level in tree
	 * @param b_i interference reduced channel value
	 * @param row_offset offset value to read from inPreCalculatedr_s
	 */
	inline void get_b_i(unsigned int *p, int d, complex<T> *b_i);
	/**
	 * \brief gets the current sphere radius
	 *
	 * @param xML current ML path
	 * @param lambdaNonML metrics of hypothesis
	 * @param path path in tree
	 * @param tree_level level in tree
	 * @param lambdaML best metric found so far
	 * @return current sphere radius
	 */
	inline T determine_radius(bool *xML, T *lambdaNonML, unsigned int *path, unsigned int tree_level, T lambdaML);
	/**
	 * \brief updates xML, lambdaML, lambdaNonML when a leaf node is reached
	 *
	 * @param new_point path to the leaf node
	 * @param new_value metric of the leaf node
	 * @param xML current ML path
	 * @param lambdaML current ML metric
	 * @param lambdaNonML metrics of hypothesis
	 */
	inline void update_minima(unsigned int *new_point, T & new_value, bool *xML, T *lambdaML, T *lambdaNonML);
	/**
	 * \brief initialisation of min0, min1 before tree search
	 *
	 * @param min0 minima for bits 0
	 * @param min1 minima for bits 1
	 */
	inline void init_minima(T *lambdaNonML, T &lambdaML);
	/**
	 * \brief computes AppPerLevel
	 *
	 * @param AppLlr apriori LLRs of current transmission vector
	 */
	inline void getAppPerLevel(T *AppLlr);
	/**
	 * \brief compute tree_app_llrs at the end of tree search
	 *
	 * @param xML ML path
	 * @param lambdaML ML metric
	 * @param lambdaNonML metrics of hypothesis
	 */
	inline void CalcAppLlr(bool *xML, T lambdaML, T *lambdaNonML);

	/**
	 * \brief clip outExtrinsicLlr to LlrClipValue
	 */
//	inline void ClipLlrs(Vec<T> & LlrVec);

	/**
	 * \brief subtracts the maximum value of the table SigmaSymbolTable from each element.
	 *
	 * This function assures that every element of SigmaSymbolTable is smaller or equal zero.
	 * Therefore the monotonous increase of the tree search metrics is not disturbed by the bias reduction.
	 */
//	inline void MonotonizeSigmaSymbolTable();

};

#endif /* MIMO_DETECTION_SPHERE_H_ */


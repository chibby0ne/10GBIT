//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * mimo_detection_sphere_param.h
 *
 *  Created on: 03.02.2011
 *      Author: gimmler
 */

#ifndef MIMO_DETECTION_SPHERE_PARAM_H_
#define MIMO_DETECTION_SPHERE_PARAM_H_

#include "../modules/base/base_param.h"
#include "../modules/transmit/mapper_share.h"

using namespace cse_lib;

///parameter class for mimo sphere decoder
class Mimo_Detection_Sphere_Parameter: public Base_Parameter {
public:
	Mimo_Detection_Sphere_Parameter() {
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Mimo_Detection_Sphere_Parameter() {
	}

	///unique name: Mimo_Detection_Sphere
	static string Unique_ID() {
		return "Mimo_Detection_Sphere";
	}

	///number of transmit/receive antennas
	Param<unsigned int> num_antennas;
	///modulation type
	Param<Mapper_Share::MAPPINGS_ENUM> modulation;
	///noise variance
	Param<double> sigma_square;
	///criteria for QR decomposition
	enum PREPROCESSING {
		ZERO_FORCING = 0, MMSE
	};
	///preprocessing optimization criterion
	Param<PREPROCESSING> preproc_criterion;
	///have the diagonal elements been normalized in the QR decomposition?
	Param<bool> diagonale_normalization;
	///number of outer iterations (over channel decoder)
	Param<int> num_outer_iterations;
	///number of inner iterations the channel decoder runs in each outer iteration
	Param<int> num_inner_iterations;

	///sphere radius
	Param<double> sphere_radius;
	///extrinsic scaling factor
	Param<double> extrinsic_scaling_factor;
	///clipping value for output LLRs
	Param<double> output_clipping_value;

	/**
	 * \brief apriori feed in for SphereDecoder
	 *
	 * In a soft-in sphere decoder the apriori information has to be incorporated into the metric calculation:
	 * \f[ \lambda = \|y - Hs\|^2 + \sum log P(x_{ij}) \f]
	 * There are several ways how to obtain the (approximative) log probabilities.
	 */
	enum APRIORI_FEED_IN {
		///straight forward implementation \f$ \frac{1}{2} x_{ij} L^a(x_{ij}) \f$
		NORMAL = 0,
		///add the full information if \f$ x_{ij} \f$ corresponds to 1
		ONLY_ONES,
		///ensures a positive increase \f$ \|L^a(x_{ij}\| + \frac{1}{2} x_{ij} L^a(x_{ij}) \f$
		POSITIVE,
		///computes the actual log probabilities out of the LLR values, ATTENTION: high complexity and NOT TESTED!!!!!!!!!!
		OPTIMAL
	};
	///apriori feed in, determines how tree_apriori_per_symbol is computed
	Param<APRIORI_FEED_IN> apriori_feedin;

	///sphere constraint for sphere decoder
	enum SPHERE_CONSTRAINT {
		///sphere constraint stays constant during tree search and is set to the sphere radius
		CONSTANT_SPHERE = 0,
		///sphere constraint is initialized with sphere radius, at each point in the tree the radius is determined by the minima found so far (see ETH paper)
		SPHERE_SHRINKING,
		///sphere constraint is the maximum value of the found minima
		EASY_SHRINKING,
		///a tree level dependent sphere constraint which is determined by the apriori information \f$ AppRadiusOffset - AppScalingFactor * \sum \| L^a(x_{ij}\|)\f$
		AACS,
		///like SPHERE_SHRINKING, however, the distance between the ML solution and the hypothesis cannot be larger than LlrClipValue
		LLR_CLIPPING,
		///constraint is determined by the maximum value in a list of K best solutions found so far
		LIST_BASED
	};
	///sphere constraint
	Param<SPHERE_CONSTRAINT> sphere_constraint;

	///tree search strategy
	enum TREE_SEARCH_STRATEGY {
		///real depth-first, requires enumeration
		DEPTH_FIRST = 0,
		///all children of a node are computed, no enumeration. this is faster than full sorting.
		FAST,
		SIC_HO,
		SIC_SO
	};
	///search strategy
	Param<TREE_SEARCH_STRATEGY> search_strategy;
	///enum for different enumeration techniques in tree search
	enum ENUMERATION_TECHNIQUE {
		NONE = 0, FULL_SORTING, CHANNEL_BASED, APRIORI_BASED, HYBRID
	};
	///enumeration technique
	Param<ENUMERATION_TECHNIQUE> enumeration_technique;
	///Threshold to change between CHANNEL_BASED and APRIORI_BASED in the HYBRID scheme
	Param<double> hybrid_threshold;

	//parameters for SPHERE_CONSTRAINT APP_LEVEL
	///parameter for SPHERE_CONSTRAINT APP_LEVEL
	Param<double> AppRadiusOffset;
	///parameter for SPHERE_CONSTRAINT APP_LEVEL
	Param<double> AppScalingFactor;

	//parameters for sphere constraint LIST_BASED
	///list size for LIST_BASED sphere constraint
	Param<unsigned int> list_size;

protected:
	///set all parameters to default values
	/**
	 * Default values:
	 *  - num_antennas				: 4
	 *  - sigma_square				: 0.1
	 *  - sphere_radius				: 2.0
	 *  - extrinsic_scaling_factor	: 1.0
	 *  - output_clipping_value		: 2.0
	 *  - diagonale_normalization	: false
	 *  - num_outer_iterations		: 1
	 *  - num_inner_iterations		: 2
	 *  - modulation				: MAP_16_QAM
	 *  - preproc_criterion			: MMSE
	 *  - apriori_feedin			: NORMAL
	 *  - sphere_constraint			: CONSTANT_SPHERE
	 *  - AppRadiusOffset			: 1.0
	 *  - AppScalingFactor			: 1.0
	 *  - list_size					: 1
	 *  - search_strategy			: FAST
	 *  - enumeration_technique		: NONE
	 *  - hybrid_threshold			: 1.0
	 */
	virtual void Set_Default_Values() {
		num_antennas.Init(4, "num_antennas", param_list_);
		sigma_square.Init(0.1, "sigma_square", param_list_);
		sphere_radius.Init(2.0, "sphere_radius", param_list_);
		extrinsic_scaling_factor.Init(1.0, "extrinsic_scaling_factor", param_list_);
		output_clipping_value.Init(2.0, "output_clipping_value", param_list_);
		diagonale_normalization.Init(false, "diagonale_normalization", param_list_);
		num_outer_iterations.Init(1, "num_outer_iterations", param_list_);
		num_inner_iterations.Init(2, "num_inner_iterations", param_list_);

		modulation.Init(Mapper_Share::MAP_16_QAM, "modulation", param_list_);
		modulation.Link_Value_String(Mapper_Share::MAP_BPSK, "MAP_BPSK");
		modulation.Link_Value_String(Mapper_Share::MAP_QPSK, "MAP_QPSK");
		modulation.Link_Value_String(Mapper_Share::MAP_QPSK_PI4, "MAP_QPSK_PI4");
		modulation.Link_Value_String(Mapper_Share::MAP_8_PSK, "MAP_8_PSK");
		modulation.Link_Value_String(Mapper_Share::MAP_16_APSK, "MAP_16_APSK");
		modulation.Link_Value_String(Mapper_Share::MAP_32_APSK, "MAP_32_APSK");
		modulation.Link_Value_String(Mapper_Share::MAP_16_QAM, "MAP_16_QAM");
		modulation.Link_Value_String(Mapper_Share::MAP_64_QAM, "MAP_64_QAM");
		modulation.Link_Value_String(Mapper_Share::MAP_256_QAM, "MAP_256_QAM");
		modulation.Link_Value_String(Mapper_Share::MAP_1024_QAM, "MAP_1024_QAM");
		modulation.Link_Value_String(Mapper_Share::MAP_4096_QAM, "MAP_4096_QAM");
		modulation.Link_Value_String(Mapper_Share::MAP_16_QAM_ITPP, "MAP_16_QAM_ITPP");
		modulation.Link_Value_String(Mapper_Share::MAP_64_QAM_ITPP, "MAP_64_QAM_ITPP");

		preproc_criterion.Init(MMSE, "preprocessing_criterion", param_list_);
		preproc_criterion.Link_Value_String(ZERO_FORCING, "ZERO_FORCING");
		preproc_criterion.Link_Value_String(MMSE, "MMSE");

		apriori_feedin.Init(NORMAL, "apriori_feedin", param_list_);
		apriori_feedin.Link_Value_String(NORMAL, "NORMAL");
		apriori_feedin.Link_Value_String(ONLY_ONES, "ONLY_ONES");
		apriori_feedin.Link_Value_String(POSITIVE, "POSITIVE");
		apriori_feedin.Link_Value_String(OPTIMAL, "OPTIMAL");

		sphere_constraint.Init(CONSTANT_SPHERE, "sphere_constraint", param_list_);
		sphere_constraint.Link_Value_String(CONSTANT_SPHERE, "CONSTANT_SPHERE");
		sphere_constraint.Link_Value_String(SPHERE_SHRINKING, "SPHERE_SHRINKING");
		sphere_constraint.Link_Value_String(EASY_SHRINKING, "EASY_SHRINKING");
		sphere_constraint.Link_Value_String(AACS, "AACS");
		sphere_constraint.Link_Value_String(LLR_CLIPPING, "LLR_CLIPPING");
		sphere_constraint.Link_Value_String(LIST_BASED, "LIST_BASED");

		//		AppOffsetIteration.Init(1, "AppOffsetIteration", param_list_);
		AppRadiusOffset.Init(1.0, "AppRadiusOffset", param_list_);
		AppScalingFactor. Init(1.0, "AppScalingFactor", param_list_);
		list_size.Init(1, "list_size", param_list_);

		search_strategy.Init(FAST, "search_strategy", param_list_);
		search_strategy.Link_Value_String(DEPTH_FIRST, "DEPTH_FIRST");
		search_strategy.Link_Value_String(FAST, "FAST");
		search_strategy.Link_Value_String(SIC_HO, "SIC_HO");
		search_strategy.Link_Value_String(SIC_SO, "SIC_SO");

		enumeration_technique.Init(NONE, "enumeration_technique", param_list_);
		enumeration_technique.Link_Value_String(NONE, "NONE");
		enumeration_technique.Link_Value_String(FULL_SORTING, "FULL_SORTING");
		enumeration_technique.Link_Value_String(CHANNEL_BASED, "CHANNEL_BASED");
		enumeration_technique.Link_Value_String(APRIORI_BASED, "APRIORI_BASED");
		enumeration_technique.Link_Value_String(HYBRID, "HYBRID");

		hybrid_threshold.Init(1.0, "hybrid_threshold", param_list_);

	}

};

#endif /* MIMO_DETECTION_SPHERE_PARAM_H_ */


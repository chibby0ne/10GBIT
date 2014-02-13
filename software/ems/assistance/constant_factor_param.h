//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * interleaving_table_based.h
 *
 *  Created on: 15.02.2011
 *      Author: gimmler
 */

#ifndef CONSTANT_FACTOR_PARAM_H_
#define CONSTANT_FACTOR_PARAM_H_

#include "../modules/base/base_param.h"

using namespace cse_lib;

///parameter class for table based interleaving
class Constant_Factor_Parameter: public Base_Parameter {
public:
	Constant_Factor_Parameter() {
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Constant_Factor_Parameter() {
	}

	///unique name: Interleaving_Table_Based
	static string Unique_ID() {
		return "Constant_Factor";
	}

	Param< float > factor;

protected:
	///set all parameters to default values
	virtual void Set_Default_Values() {
		factor.Init(1, "factor", param_list_);
	}

};

#endif /* CONSTANT_FACTOR_PARAM_H_ */

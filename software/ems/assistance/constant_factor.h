//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/// \file
/// \brief  Definition of the constant factor module.
/// \author Christina Gimmler
/// \date   09.01.2012
//

#ifndef CONSTANT_FACTOR_H_
#define CONSTANT_FACTOR_H_

#include <complex>
#include "constant_factor_iface.h"
#include "constant_factor_param.h"

/// Multiply every element of a buffer with a constant factor
template<class T>
class Constant_Factor: public Constant_Factor_Interface<T> , public Constant_Factor_Parameter {

public:

	Constant_Factor() {
		buffer_length = 0;
	}
	virtual ~Constant_Factor() {
	}

	int Run();

protected:
	///initialize all values
	void Init();
	///resize all memories
	void Alloc_Mem();

private:

	int buffer_length;

};

#endif // CONSTANT_FACTOR_H_

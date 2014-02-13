//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  SISO mapper class definition.
/// \author Matthias Alles
/// \date   2010/12/10
//

#ifndef MAPPER_H_
#define MAPPER_H_

#include <complex>
#include "mapper_iface.h"
#include "mapper_param.h"
#include "mapper_share.h"


namespace cse_lib {

/// Flexible SISO Mapper, mapping bits to complex modulation symbols.
/**
 * The Mapper supports BPSK, QPSK, 8PSK, 16-APSK, 32-APSK, 16-QAM, 64-QAM, 256-QAM, etc.
 *
 * \ingroup modules
 */
class Mapper: public Mapper_Interface,
              public Mapper_Parameter,
              public Mapper_Share
{

public:

	Mapper();
	virtual ~Mapper();

	/// Process one block
	int Run();

private:

	void Alloc_Mem();
	void Init();

	unsigned int bits_per_symbol_;
};
}
#endif // MAPPER_H_

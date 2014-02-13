//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * mimo_channel_rayleigh.h
 *
 *  Created on: 03.02.2011
 *      Author: gimmler
 */

#ifndef MIMO_CHANNEL_RAYLEIGH_H_
#define MIMO_CHANNEL_RAYLEIGH_H_

#include "itpp/itbase.h"

#include "mimo_channel_rayleigh_iface.h"
#include "mimo_channel_rayleigh_param.h"

using namespace itpp;
using namespace cse_lib;

/**
 * \brief Rayleigh fading mimo channel
 *
 * This class realizes a Rayleigh fading MIMO channel with a symmetric antenna constellation.
 * The channel fades independently after an adjustable number of channel uses.
 *
 * The Run function transmits the modulation symbols of the input buffer over a Rayleigh fading MIMO channel.
 *
 * Channel Model: \f[y = \sqrt{\frac{1}{M_T}} H s + \sqrt{\frac{1}{E_S/N_0}} n\f]
 * The received channel values and all channel realizations are written to the according output buffers.
 */
template<class T>
class Mimo_Channel_Rayleigh: public Mimo_Channel_Rayleigh_Interface<T> , public Mimo_Channel_Rayleigh_Parameter {
public:
	Mimo_Channel_Rayleigh();
	virtual ~Mimo_Channel_Rayleigh();
	///process one block of input data
	int Run();

protected:
	///compute derived parameters and resize memories
	void Init();
	///resize all memories
	void Alloc_Mem();

private:

	///scaling factor for Gaussian noise
	T noise_factor;
	///scaling factor for channel matrix
	T channel_factor;
	///number of channel uses per data block
	int num_channel_uses_per_block;
	///number of modulation symbols per data block
	int num_symbols_per_block;

};

#endif /* MIMO_CHANNEL_RAYLEIGH_H_ */


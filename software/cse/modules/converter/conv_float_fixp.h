//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the float to fixpoint-converter module.
/// \author Matthias Alles
/// \date   2010/12/13
//

#ifndef CONV_FLOAT_FIXP_H_
#define CONV_FLOAT_FIXP_H_

#include <complex>
#include "conv_float_fixp_iface.h"
#include "conv_float_fixp_param.h"


namespace cse_lib {


/// Convert a Buffer of type float or complex<float> to fixpoint types with saturation.
/**
 * Convert the input Buffer of type float or complex<float> to fixpoint values of type
 * int or complex<int>. The values are scaled and saturated to two's complement value
 * range.
 *
 * \ingroup modules
 */
template<class T_in, class T_out>
class Converter_Float_Fixpoint: public Converter_Float_Fixpoint_Interface<T_in, T_out>,
                                public Converter_Float_Fixpoint_Parameter
{

public:

	Converter_Float_Fixpoint() { }
	virtual ~Converter_Float_Fixpoint() { }

	int Run();

private:

	/// Maximum magnitude representable with bw_output (in two's complement arithmetic).
	int max_magnitude_;

	/// Scaling factor required because of fractional bits (nun_bits_fract).
	unsigned int mult_factor_;

	/// Convert and saturate the float values to fixpoint values.
	void Convert(Buffer<complex<float> > &input,
	             Buffer<complex<int> >   &output);

	/// Convert and saturate the float values to fixpoint values.
	void Convert(Buffer<float> &input,
	             Buffer<int>   &output);

};
}
#endif // CONV_FLOAT_FIXP_H_

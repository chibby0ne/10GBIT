//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  SISO mapper shared functionality class definition.
/// \author Matthias Alles
/// \date   2010/12/10
//

#ifndef MAPPER_SHARE_H_
#define MAPPER_SHARE_H_

#include <complex>
#include "../../assistance/buffer.h"
#include "../base/cse_namespace.h"

#include "itpp/itbase.h"

using namespace itpp;

namespace cse_lib {


/// Functions that can be shared between multiple mappers.
/**
 * In particular this module contains the mappings for many modulation schemes
 * in form of look-up-tables (LUTs).
 * The LUTs are of such a form, that the index address corresponds to the bit
 * sequence of the contained constellation point at that address. <b>The MSB of
 * the LUT address corresponds to the bit that is being transmitted first!</b>
 *
 * \ingroup share
 */
class Mapper_Share
{

public:

	/// Available modulation schemes.
	enum MAPPINGS_ENUM
	{
		MAP_BPSK,      ///< Standard BPSK mapping (1 -> -1, 0 -> 1).
		MAP_QPSK,      ///< Standard QPSK mapping, compliant with, e.g., DVB-S2.
		MAP_QPSK_PI4,  ///< QPSK mapping rotated by Pi/4 (i.e., constellation points on axes),
		               ///< compliant with GMR Release 2 and 3.
		MAP_8_PSK,     ///< 8-PSK mapping, compliant with DVB-S2 or GMR Release 2 and 3.
		MAP_16_APSK,   ///< 16-APSK mapping, compliant with DVB-S2 or GMR Release 2 and 3.
		MAP_32_APSK,   ///< 32-APSK mapping, compliant with DVB-S2 or GMR Release 2 and 3.
		MAP_16_QAM,    ///< 16-QAM gray mapping.
		MAP_64_QAM,    ///< 64-QAM gray mapping.
		MAP_256_QAM,   ///< 246-QAM gray mapping.
		MAP_1024_QAM,  ///< 1024-QAM gray mapping.
		MAP_4096_QAM,  ///< 4096-QAM gray mapping.
		MAP_16_QAM_ITPP, ///< 16-QAM gray mapping ITPP style.
		MAP_64_QAM_ITPP  ///< 64-QAM gray mapping ITPP style.
	};

	/// Function which holds the corresponding strings to the enum type
	static DS_Enum_String<MAPPINGS_ENUM> Mapping_Enum_Str()
	{
		DS_Enum_String<MAPPINGS_ENUM> data_struct;
		data_struct(MAP_BPSK    , "MAP_BPSK"    )
		           (MAP_QPSK    , "MAP_QPSK"    )
		           (MAP_QPSK_PI4, "MAP_QPSK_PI4")
		           (MAP_8_PSK   , "MAP_8_PSK"   )
		           (MAP_16_APSK , "MAP_16_APSK" )
		           (MAP_32_APSK , "MAP_32_APSK" )
		           (MAP_16_QAM  , "MAP_16_QAM"  )
		           (MAP_64_QAM  , "MAP_64_QAM"  )
		           (MAP_256_QAM , "MAP_256_QAM" )
		           (MAP_1024_QAM, "MAP_1024_QAM")
		           (MAP_4096_QAM, "MAP_4096_QAM")
		           (MAP_16_QAM_ITPP  , "MAP_16_QAM_ITPP"  )
		           (MAP_64_QAM_ITPP  , "MAP_64_QAM_ITPP"  );

		return data_struct;
	}

protected:

	/// Dump a message about an unsupported constellation
	/**
	 * This is in particular useful, if demappers are not fully generic and
	 * support only a subset of the constellations as defined in this class.
	 */
	void Dump_Wrong_Constellation_Message(string instance_name, string mapping);

	/// Translate modulation/mapping scheme into bits_per_symbol.
	unsigned int Get_Bits_Per_Symbol(enum MAPPINGS_ENUM mapping);

	/// Set the map_lut_ content according to the chosen parameters.
	/**
	 * \return Return the power normalization factor.  This factor is <= 1 and
	 * was used to normalize the symbols within the mapping LUT to energy 1 on
	 * average.
	 */
	float Get_Mapping_LUT(enum MAPPINGS_ENUM mapping,
	                     float psk_radius_ratio_r2,
	                     float psk_radius_ratio_r3);

	/// Normalize the overall power of all symbols to 1.
	float Normalize_Power();

	/// Perform mapping of input bits to symbols with help of the mapping LUT (map_lut_).
	/**
	 * Padding is not supported that means the number of bits_per_symbol and the number of input_bits
	 * have to be an integral multiple.
	 */
	void Bit_Mapping(unsigned int             bits_per_symbol,
	                 Buffer<unsigned int>    &input_bits,
	                 Buffer<complex<float> > &output_symb);

	/// Generic function for QAM modulation
	void Get_Mapping_LUT_xQAM(unsigned int bits_per_symbol);

	/// Normal BPSK modulation
	void Get_Mapping_LUT_BPSK();

	/// Normal QPSK modulation
	void Get_Mapping_LUT_QPSK();

	/// Pi/4-QPSK modulation
	void Get_Mapping_LUT_QPSK_PI4();

	/// 8-PSK modulation according to DVB-S2 standard.
	void Get_Mapping_LUT_8PSK();

	/// 16-APSK modulation with constellation points according to DVB-S2 standard.
	void Get_Mapping_LUT_16APSK(float psk_radius_ratio_r2);

	/// 32-APSK modulation with constellation points according to DVB-S2 standard.
	void Get_Mapping_LUT_32APSK(float psk_radius_ratio_r2, float psk_radius_ratio_r3);

	/// LUT of constellation points for currently used mapping.
	Buffer<complex<float> > map_lut_;

	/// ITPP Matrix where row k contains the constellation points with '0' in bit position k.
	imat S0;

	/// ITPP Matrix where row k contains the constellation points with '1' in bit position k.
	imat S1;

	/// Other ITPP vars // TODO cleanup
	int k, M, L;
	bmat bitmap;
	Vec< complex< int > > symbols;
	ivec bits2symbols;

	bmat graycode(int m);
	/*int  bin2dec(const bvec &inbvec, bool msb_first);*/
	void calculate_softbit_matrices();
	void Get_Mapping_LUT_ITPP_xQAM(unsigned int bits_per_symbol);


};
}
#endif // MAPPER_SHARE_H_

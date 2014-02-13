/*
 * WPAN_chain.cpp
 *
 *  Created on: 09.11.2012
 *      Author: schlaefer
 */

#include "../cse/include/cse_lib.h"
//#include "ldpc_enc/enc_ldpc_ieee_802_11ad.h"
#include "ldpc_dec/dec_ldpc_ieee_802_3an.h"

#include <iostream>
#include <fstream>

using namespace cse_lib;

int main(int argc, char* argv[]) {

	string configfile = "config.xml";
	if (argc >= 2)
		configfile = argv[1];

	string rfilename = "./results/";
	if (argc >= 3)
		rfilename = rfilename + argv[2];
	else
		rfilename = rfilename + "results.xml";

	Manage_Module_Config xml_config(configfile);
	Manage_Module_Result xml_result(rfilename, "WPAN");

	// Instantiate modules
	Source_Bits source_bits;

//	Encoder_LDPC_IEEE_802_11ad encoder;

	Mapper mapper;

	Converter_Float_Fixpoint<float, int> converter;

	Channel_AWGN channel;

	Demapper demapper;

    Decoder_LDPC_IEEE_802_3an<int> decoder;

	Statistics_Error_Rates<2> error_rates_decoding;
	error_rates_decoding.instance_name("error_rates_decoding");


	// Connect modules
//	encoder.input_bits(source_bits.output_bits());
//	mapper.input_bits(encoder.output_bits());
	mapper.input_bits(source_bits.output_bits());
	channel.input_symb(mapper.output_symb());
	demapper.input_symb(channel.output_symb());
	converter.input(demapper.output_bits_llr());

	decoder.input_bits_llr(converter.output());

	error_rates_decoding.input_bits_ref(source_bits.output_bits());
	error_rates_decoding.input_bits(decoder.output_bits());

	int result;
	do {
		// Configure modules
		xml_config.Configure_Module(source_bits);
//		xml_config.Configure_Module(encoder);
		xml_config.Configure_Module(mapper);
		xml_config.Configure_Module(channel);
		xml_config.Configure_Module(demapper);
		xml_config.Configure_Module(converter);
		xml_config.Configure_Module(decoder);
		xml_config.Configure_Module(error_rates_decoding);

		//RNG_reset();

		do {
			source_bits.Run();
//			encoder.Run();
			mapper.Run();

			channel.Run();

			demapper.Run();
			converter.Run();
			decoder.Run();


			result = error_rates_decoding.Run(); 
		} while (result == 0);

		xml_result.Create_Iteration_Value_Result_Point(xml_config);  // Create a new iteration value XML tree to store the modules results (uses xml_config to as a template for results)
		xml_result.Insert_Results_From_Module(error_rates_decoding); // Insert the results from a Module into the current working tree
		xml_result.Insert_Results_From_Module(decoder); // Insert the results from a Module into the current working tree
		xml_result.Write_Current_State(); // Write the current XML result tree into the current working tree
	} while (xml_config.Update_To_Next_Iter() == 0); // Update configuration instance with the next iteration point

	return 0;
}


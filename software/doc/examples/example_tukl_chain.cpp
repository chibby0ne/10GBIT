//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  CSE example chain.
/// \author Matthias Alles
/// \date   2011/01/18
//


#include "cse_lib.h"

using namespace cse_lib;

int main(int argc, char *argv[])
{

	Manage_Module_Config xml_config("config.xml");
	Manage_Module_Result xml_result("result.xml",argv[0]);

	// Instantiate modules
	Source_Bits source_bits;
	Encoder_LDPC_GMR enc_ldpc;
	Interleaving_GMR<unsigned int> interleaver;
	interleaver.instance_name("Interleaver_GMR");
	Mapper mapper;
	Channel_AWGN channel_awgn;
	Demapper demapper;
	Converter_Float_Fixpoint<float, int> converter;
	Interleaving_GMR<int> deinterleaver;
	deinterleaver.instance_name("Deinterleaver_GMR");
	Decoder_LDPC_GMR dec_ldpc;
	Statistics_Error_Rates<2> error_rates;



	// Connect modules
	enc_ldpc.input_bits(source_bits.output_bits());
	interleaver.input(enc_ldpc.output_bits());
	mapper.input_bits(interleaver.output());
	channel_awgn.input_symb(mapper.output_symb());
	demapper.input_symb(channel_awgn.output_symb());
	converter.input(demapper.output_bits_llr());
	deinterleaver.input(converter.output());
	dec_ldpc.input_bits_llr(deinterleaver.output());
	error_rates.input_bits(dec_ldpc.output_bits());
	error_rates.input_bits_ref(source_bits.output_bits());
 
	int result;

	do
	{
		// Configure modules
		xml_config.Configure_Module(source_bits);
		xml_config.Configure_Module(enc_ldpc);
		xml_config.Configure_Module(interleaver);
		xml_config.Configure_Module(mapper);
		xml_config.Configure_Module(channel_awgn);
		xml_config.Configure_Module(converter);
		xml_config.Configure_Module(demapper);
		xml_config.Configure_Module(deinterleaver);
		xml_config.Configure_Module(dec_ldpc);
		xml_config.Configure_Module(error_rates);

		do
		{
			source_bits.Run();
			enc_ldpc.Run();
			interleaver.Run();
			mapper.Run();

			channel_awgn.Run();

			demapper.Run();
			converter.Run();
			deinterleaver.Run();
			dec_ldpc.Run();

			result = error_rates.Run();
		}
		while(result == 0);

		xml_result.Create_Iteration_Value_Result_Point(xml_config);
		xml_result.Insert_Results_From_Module(demapper);
		xml_result.Insert_Results_From_Module(error_rates);
		xml_result.Write_Current_State();
	}
	while(xml_config.Update_To_Next_Iter() == 0);

	return 0;
}

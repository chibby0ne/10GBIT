//
//  Copyright (C) 2011 TU Kaiserslautern, Microelectronic Systems Design Research Group
//
/*
 * interleaving_table_based.h
 *
 *  Created on: 15.02.2011
 *      Author: gimmler
 */

#ifndef INTERLEAVING_TABLE_BASED_PARAM_H_
#define INTERLEAVING_TABLE_BASED_PARAM_H_

using namespace cse_lib;

///parameter class for table based interleaving
class Interleaving_Table_Based_Parameter: public Base_Parameter {
public:
	Interleaving_Table_Based_Parameter() {
		instance_name(Unique_ID());
		Set_Default_Values();
	}

	virtual ~Interleaving_Table_Based_Parameter() {
	}

	///unique name: Interleaving_Table_Based
	static string Unique_ID() {
		return "Interleaving_Table_Based";
	}

	///filename of stored interleaver table
	Param<string> file_name;

	///choice of available interleaver tables
	enum INTERLEAVER_TABLE {
		NO_TABLE,///no given table, use given file_name
		LTE_2304,///LTE compliant interleaver of length 2304
		RP_2304
	///Random Prime Interleaver of length 2304
	};
	///chosen interleaver table
	Param<INTERLEAVER_TABLE> interleaver_table;
	///path to given interleaver tables
	Param<string> table_path;

	///modi for interleavers
	enum INTERLEAVING_MODE {
		INTERLEAVE, ///Perform Interleaving when calling Run().
		DEINTERLEAVE
	///Perform Deinterleaving when calling Run().
	};

	/// Use module as interleaver or deinterleaver?
	Param<INTERLEAVING_MODE> mode;

protected:
	///set all parameters to default values
	/**
	 * Default values:
	 *  - file_name			: interleaver_table.txt
	 *  - mode				: INTERLEAVE
	 *  - interleaver_table	: NO_TABLE
	 *  - table_path		: /net/user/r1/unix/gimmler/projekte/emssim/library/interleaving/interleaver_tables/
	 */
	virtual void Set_Default_Values() {
		file_name.Init("interleaver_table.txt", "file_name", param_list_);
		mode.Init(INTERLEAVE, "mode", param_list_);
		mode.Link_Value_String(INTERLEAVE, "INTERLEAVE");
		mode.Link_Value_String(DEINTERLEAVE, "DEINTERLEAVE");
		interleaver_table.Init(NO_TABLE, "interleaver_table", param_list_);
		interleaver_table.Link_Value_String(NO_TABLE, "NO_TABLE");
		interleaver_table.Link_Value_String(LTE_2304, "LTE_2304");
		interleaver_table.Link_Value_String(RP_2304, "RP_2304");
		//TODO maybe a relative path considering projects structure would be better?
		table_path.Init("/net/user/r1/unix/gimmler/projekte/emssim/library/interleaving/interleaver_tables/", "table_path", param_list_);
//		table_path.Init("../../library/interleaving/interleaver_tables/", "table_path", param_list_);
	}

};

#endif /* INTERLEAVING_TABLE_BASED_PARAM_H_ */

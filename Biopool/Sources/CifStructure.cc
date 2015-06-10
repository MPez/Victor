/* 
 * File:   CifStructure.cc
 * Author: marco
 * 
 * Created on 1 giugno 2015, 11.36
 */

#include <iostream>
#include <algorithm>

#include <IoTools.h>

#include "CifStructure.h"

using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

CifStructure::CifStructure(istream& input, ostream& output) :
input(input), output(output) {
    header = "_entry.id";
    model = "pdbx_PDB_model_num ";
    helix = "_struct_conf.";
    helixStart = "beg_auth_seq_id ";
    helixEnd = "end_auth_seq_id ";
    helixChainId = "beg_auth_asym_id ";
    atom = "_atom_site.";
    residueNum = "auth_seq_id ";
    atomId = "id ";
    residueIns = "pdbx_PDB_ins_code ";
    tempFactor = "B_iso_or_equiv ";
    atomName = "auth_atom_id ";
    residueName = "auth_comp_id ";
    x = "Cartn_x ";
    y = "Cartn_y ";
    z = "Cartn_z ";
    chain = "auth_asym_id ";
    sheet = "_struct_sheet.";
    sheetOrder = "_struct_sheet_order.";
    sheetRange = "_struct_sheet_range.";
    sheetHbond = "_pdbx_struct_sheet_hbond.";
    sheetStart = "beg_auth_seq_id ";
    sheetEnd = "end_auth_seq_id ";
    sheetChainId = "beg_auth_asym_id ";

    atomGroupParsed = false;
    helixGroupParsed = false;
    sheetGroupParsed = false;
    sheetOrderGroupParsed = false;
    sheetRangeGroupParsed = false;
    sheetHboundgroupParsed = false;
}

CifStructure::~CifStructure() {
}

/**
 * Returns the correct collection by group name
 * @param name name of the CIF group
 * @return reference to the collection
 */
vector<string>& CifStructure::getGroup(string name) {
    if (name == "atom") {
	return atomGroup;
    } else if (name == "helix") {
	return helixGroup;
    } else if (name == "sheet") {
	return sheetGroup;
    } else if (name == "sheet order") {
	return sheetOrderGroup;
    } else if (name == "sheet range") {
	return sheetRangeGroup;
    } else if (name == "sheet hbond") {
	return sheetHbondGroup;
    }
}

/**
 * Returns the tag by name
 * @param name name of tag
 * @return CIF tag
 */
string CifStructure::getTag(string name) {
    if (name == "header") {
	return header;
    } else if (name == "atom") {
	return atom;
    } else if (name == "residue num") {
	return residueNum;
    } else if (name == "atom id") {
	return atomId;
    } else if (name == "residue ins") {
	return residueIns;
    } else if (name == "x") {
	return x;
    } else if (name == "y") {
	return y;
    } else if (name == "z") {
	return z;
    } else if (name == "bfac") {
	return tempFactor;
    } else if (name == "atom name") {
	return atomName;
    } else if (name == "residue name") {
	return residueName;
    } else if (name == "helix") {
	return helix;
    } else if (name == "helix start") {
	return helixStart;
    } else if (name == "helix end") {
	return helixEnd;
    } else if (name == "helix chain") {
	return helixChainId;
    } else if (name == "model") {
	return model;
    } else if (name == "chain") {
	return chain;
    } else if (name == "sheet") {
	return sheet;
    } else if (name == "sheet order") {
	return sheetOrder;
    } else if (name == "sheet range") {
	return sheetRange;
    } else if (name == "sheet hbond") {
	return sheetHbond;
    } else if (name == "sheet start") {
	return sheetStart;
    } else if (name == "sheet end") {
	return sheetEnd;
    } else if (name == "sheet chain") {
	return sheetChainId;
    }
}

/**
 * Returns the column number of the field
 * @param name name of the group
 * @param field name of the field
 * @return field column number
 */
int CifStructure::getGroupColumnNumber(string name, string field) {
    int col = -1;
    vector<string>& group = getGroup(name);
    string tag = getTag(field);
    for (vector<string>::iterator it = group.begin(); it != group.end(); it++) {
	if (*it == tag) {
	    col = it - group.begin();
	    break;
	}
    }
    return col;
}

/**
 * Returns the field of the line at the columnNum column
 * @param name name of the group
 * @param line line of the CIF
 * @param columnNum number of column
 * @return field at columnNum column
 */
string CifStructure::getGroupField(string name, string& line, int columnNum) {
    istringstream iss(line);
    vector<string>& group = getGroup(name);
    vector<string> fields;
    string field;
    for (unsigned int i = 0; i < group.size(); i++) {
	iss >> field;
	fields.push_back(field);
    }
    if (columnNum != -1) {
	return fields[columnNum];
    } else {
	return "?";
    }
}

/**
* Returns the field present in the line
* @param line line of the CIF
* @return field of the line
*/
string CifStructure::getInlineField(string& line) {
    istringstream iss(line);
    string tag, field;
    iss >> tag >> field;
    return field;
}

/**
 * Parses group of CIF fields and creates a vector with columns positions
 * @param name name of the group 
 */
void CifStructure::parseGroup(string name, string& line) {
    bool found = false;
    vector<string>& group = getGroup(name);

    // exit the function if the group name is already parsed
    if (!isGroupParsed(name)) {
	while (input) {
	    string groupName(getTag(name));
	    size_t pos = line.find(groupName);
	    if (pos != string::npos) {
		group.push_back(line.substr(pos + groupName.size(),
			line.size() - groupName.size()));
		found = true;
	    } else {
		found = false;
	    }

	    // exit the loop when the research of the fields is completed
	    if (!found && group.size() > 1) {
		setParsedFlag(name);
		break;
	    }
	    
	    line = readLine(input);
	}
    }
}

/**
 * Sets flag of the parsed group
 * @param name name of the group
 */
void CifStructure::setParsedFlag(string name) {
    if (name == "atom") {
	atomGroupParsed = true;
    } else if (name == "helix") {
	helixGroupParsed = true;
    } else if (name == "sheet") {
	sheetGroupParsed = true;
    } else if (name == "sheet hbond") {
	sheetHboundgroupParsed = true;
    } else if (name == "sheet order") {
	sheetOrderGroupParsed = true;
    } else if (name == "sheet range") {
	sheetRangeGroupParsed = true;
    }
}

/**
 * Return true if the group name is parsed, false otherwise
 * @param name name of the group
 * @return true if group is parsed, false otherwise
 */
bool CifStructure::isGroupParsed(string name) {
    if (name == "atom") {
	return atomGroupParsed;
    } else if (name == "helix") {
	return helixGroupParsed;
    } else if (name == "sheet") {
	return sheetGroupParsed;
    } else if (name == "sheet hbond") {
	return sheetHboundgroupParsed;
    } else if (name == "sheet order") {
	return sheetOrderGroupParsed;
    } else if (name == "sheet range") {
	return sheetRangeGroupParsed;
    }
}

/**
* Prints group records names into output stream.
* @param name name of the group
*/
void CifStructure::printGroup(string name) {
    vector<string>& group = getGroup(name);

    output << "loop_" << endl;

    if (isGroupParsed(name)) {
	for (vector<string>::iterator it = group.begin(); it != group.end(); it++) {
	    output << *it << endl;
	}
    } else {
	if (name == "atom") {
	    output << "_atom_site.group_PDB " << endl <<
		    "_atom_site.id " << endl <<
		    "_atom_site.type_symbol " << endl <<
		    "_atom_site.label_atom_id " << endl <<
		    "_atom_site.label_alt_id " << endl <<
		    "_atom_site.label_comp_id " << endl <<
		    "_atom_site.label_asym_id " << endl <<
		    "_atom_site.label_entity_id " << endl <<
		    "_atom_site.label_seq_id " << endl <<
		    "_atom_site.pdbx_PDB_ins_code " << endl <<
		    "_atom_site.Cartn_x " << endl <<
		    "_atom_site.Cartn_y " << endl <<
		    "_atom_site.Cartn_z " << endl <<
		    "_atom_site.occupancy " << endl <<
		    "_atom_site.B_iso_or_equiv " << endl <<
		    "_atom_site.Cartn_x_esd " << endl <<
		    "_atom_site.Cartn_y_esd " << endl <<
		    "_atom_site.Cartn_z_esd " << endl <<
		    "_atom_site.occupancy_esd " << endl <<
		    "_atom_site.B_iso_or_equiv_esd " << endl <<
		    "_atom_site.pdbx_formal_charge " << endl <<
		    "_atom_site.auth_seq_id " << endl <<
		    "_atom_site.auth_comp_id " << endl <<
		    "_atom_site.auth_asym_id " << endl <<
		    "_atom_site.auth_atom_id " << endl <<
		    "_atom_site.pdbx_PDB_model_num " << endl;
	} else if (name == "helix") {
	    output << "_struct_conf.conf_type_id " << endl <<
		    "_struct_conf.id " << endl <<
		    "_struct_conf.pdbx_PDB_helix_id " << endl <<
		    "_struct_conf.beg_label_comp_id " << endl <<
		    "_struct_conf.beg_label_asym_id " << endl <<
		    "_struct_conf.beg_label_seq_id " << endl <<
		    "_struct_conf.pdbx_beg_PDB_ins_code " << endl <<
		    "_struct_conf.end_label_comp_id " << endl <<
		    "_struct_conf.end_label_asym_id " << endl <<
		    "_struct_conf.end_label_seq_id " << endl <<
		    "_struct_conf.pdbx_end_PDB_ins_code " << endl <<
		    "_struct_conf.beg_auth_comp_id " << endl <<
		    "_struct_conf.beg_auth_asym_id " << endl <<
		    "_struct_conf.beg_auth_seq_id " << endl <<
		    "_struct_conf.end_auth_comp_id " << endl <<
		    "_struct_conf.end_auth_asym_id " << endl <<
		    "_struct_conf.end_auth_seq_id " << endl <<
		    "_struct_conf.pdbx_PDB_helix_class " << endl <<
		    "_struct_conf.details " << endl <<
		    "_struct_conf.pdbx_PDB_helix_length " << endl;
	} else if (name == "sheet") {
	    output << "_struct_sheet.id " << endl <<
		    "_struct_sheet.type " << endl <<
		    "_struct_sheet.number_strands " << endl <<
		    "_struct_sheet.details " << endl;
	} else if (name == "sheet range") {
	    output << "_struct_sheet_range.sheet_id " << endl <<
		    "_struct_sheet_range.id " << endl <<
		    "_struct_sheet_range.beg_label_comp_id " << endl <<
		    "_struct_sheet_range.beg_label_asym_id " << endl <<
		    "_struct_sheet_range.beg_label_seq_id " << endl <<
		    "_struct_sheet_range.pdbx_beg_PDB_ins_code " << endl <<
		    "_struct_sheet_range.end_label_comp_id " << endl <<
		    "_struct_sheet_range.end_label_asym_id " << endl <<
		    "_struct_sheet_range.end_label_seq_id " << endl <<
		    "_struct_sheet_range.pdbx_end_PDB_ins_code " << endl <<
		    "_struct_sheet_range.symmetry " << endl <<
		    "_struct_sheet_range.beg_auth_comp_id " << endl <<
		    "_struct_sheet_range.beg_auth_asym_id " << endl <<
		    "_struct_sheet_range.beg_auth_seq_id " << endl <<
		    "_struct_sheet_range.end_auth_comp_id " << endl <<
		    "_struct_sheet_range.end_auth_asym_id " << endl <<
		    "_struct_sheet_range.end_auth_seq_id " << endl;
	} else if (name == "sheet hbond") {
	    output << "_pdbx_struct_sheet_hbond.sheet_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_id_1 " << endl <<
		    "_pdbx_struct_sheet_hbond.range_id_2 " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_label_atom_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_label_comp_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_label_asym_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_label_seq_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_PDB_ins_code " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_auth_atom_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_auth_comp_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_auth_asym_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_1_auth_seq_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_label_atom_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_label_comp_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_label_asym_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_label_seq_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_PDB_ins_code " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_auth_atom_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_auth_comp_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_auth_asym_id " << endl <<
		    "_pdbx_struct_sheet_hbond.range_2_auth_seq_id " << endl;
	}
    }
}

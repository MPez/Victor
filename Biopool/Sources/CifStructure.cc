/* 
 * File:   CifStructure.cc
 * Author: marco
 * 
 * Created on 1 giugno 2015, 11.36
 */

#include <regex>
#include <iostream>

#include <IoTools.h>

#include "CifStructure.h"

using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

CifStructure::CifStructure(istream& input) : input(input) {
}

CifStructure::~CifStructure() {
}

/**
 * returns the correct collection by group name
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
 * returns the tag by name
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
    } else if (name == "alt id") {
        return atomAltId;
    } else if (name == "x") {
        return x;
    } else if (name == "y") {
        return y;
    } else if (name == "z") {
        return z;
    } else if (name == "bfac") {
        return tempFactor;
    } else if (name == "bfac") {
        return tempFactor;
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
 * returns the column number of the field
 * @param name name of the group
 * @param field name of the field
 * @return field column number
 */
int CifStructure::getGroupColumnNumber(string name, string field) {
    int col = -1;
    vector<string> group = getGroup(name);
    vector<string>::iterator it;
    it = find(group.begin(), group.end(), getTag(field));
    if (it != group.end()) {
        col = it - group.begin();
    }
    return col;
}

/**
 * returns the field of the line at the columnNum column
 * @param name name of the group
 * @param line line of the CIF
 * @param columnNum number of column
 * @return field at columnNum column
 */
string CifStructure::getGroupField(string name, string line, int columnNum) {
    istringstream iss(line);
    vector<string>& group = getGroup(name);
    vector<string> fields;
    string field;
    for (unsigned int i = 0; i < group.size(); i++) {
        iss >> field;
        fields.push_back(field);
    }
    return fields[columnNum];
}

/**
 * parses group of CIF fields and creates a vector with columns positions
 * @param name name of the group 
 */
void CifStructure::parseGroup(string name, string line) {
    bool found = false;
    vector<string>& group = getGroup(name);

    // exit the function if the group name is already parsed
    if (!isGroupParsed(name)) {
        while (input) {
            regex groupName(getTag(name));
            smatch match;
            if (regex_search(line, match, groupName)) {
                group.push_back(match.suffix().str());
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
    } else if (name == "sheet hbound") {
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
    } else if (name == "sheet hbound") {
        return sheetHboundgroupParsed;
    } else if (name == "sheet order") {
        return sheetOrderGroupParsed;
    } else if (name == "sheet range") {
        return sheetRangeGroupParsed;
    }
}
/*
 * File:   CifStructure.h
 * Author: marco
 *
 * Created on 1 giugno 2015, 11.36
 */

#ifndef CIFSTRUCTURE_H
#define	CIFSTRUCTURE_H


// Includes:
#include <string>
#include <vector>
#include <istream>

using std::string;
using std::istream;
using std::vector;

// Global constants, typedefs, etc. (to avoid):

namespace Victor {
    namespace Biopool {

        /**
         * Helper class used to hold information from CIF file
         */
        class CifStructure {
        public:
            /**
             * Constructor
             * @param input input file stream
             */
            CifStructure(istream& input);
            
            /**
             * Destructor
             */
            virtual ~CifStructure();
            
            /**
            * Returns the correct collection by group name.
            * @param name name of the CIF group
            * @return reference to the collection
            */
            vector<string>& getGroup(string name);
            
            /**
            * Returns the tag by name.
            * @param name name of tag
            * @return CIF tag
            */
            string getTag(string name);
            
            /**
            * Returns the column number of the field.
            * @param name name of the group
            * @param field name of the field
            * @return field column number
            */
            int getGroupColumnNumber(string name, string field);
            
            /**
            * Returns the field of the line at the columnNum column.
            * @param name name of the group
            * @param line line of the CIF
            * @param columnNum number of column
            * @return field at columnNum column
            */
            string getGroupField(string name, string line, int columnNum);
            
            /**
            * Parses group of CIF fields and creates a vector with columns positions.
            * @param name name of the group 
            */
            void parseGroup(string group, string line);
            
            /**
            * Sets flag of the parsed group
            * @param name name of the group
            */
            void setParsedFlag(string name);
            
            /**
             * Return true if the group name is parsed, false otherwise
             * @param name name of the group
             * @return true if group is parsed, false otherwise
             */
            bool isGroupParsed(string name);
            
        private:
            // CIF file
            istream& input;

            // CIF tags
            string header = "_struct_keywords.pdbx_keywords";
            string model = "pdbx_PDB_model_num";
            string helix = "_struct_conf.";
            string helixStart = "beg_auth_seq_id";
            string helixEnd = "end_auth_seq_id";
            string helixChainId = "beg_auth_asym_id";
            string atom = "_atom_site.";
            string residueNum = "auth_seq_id";
            string atomId = "id";
            string atomAltId = "label_alt_id";
            string tempFactor = "B_iso_or_equiv";
            string atomName = "auth_atom_id";
            string residueName = "auth_comp_id";
            string x = "Cartn_x";
            string y = "Cartn_y";
            string z = "Cartn_z";
            string chain = "auth_asym_id";
            string sheet = "_struct_sheet.";
            string sheetOrder = "_struct_sheet_order.";
            string sheetRange = "_struct_sheet_range.";
            string sheetHbond = "_pdbx_struct_sheet_hbond.";
            string sheetStart = "beg_auth_seq_id";
            string sheetEnd = "end_auth_seq_id";
            string sheetChainId = "beg_auth_asym_id";
            
            // collections of CIF group fields
            vector<string> atomGroup;
            vector<string> helixGroup;
            vector<string> sheetGroup;
            vector<string> sheetOrderGroup;
            vector<string> sheetRangeGroup;
            vector<string> sheetHbondGroup;

            // flags
            bool atomGroupParsed = false;
            bool helixGroupParsed = false;
            bool sheetGroupParsed = false;
            bool sheetOrderGroupParsed = false;
            bool sheetRangeGroupParsed = false;
            bool sheetHboundgroupParsed = false;
        };
    } // namespace Biopool
} // namespace Victor 

#endif	/* CIFSTRUCTURE_H */


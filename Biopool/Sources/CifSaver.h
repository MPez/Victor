/* 
 * File:   CifSaver.h
 * Author: marco
 *
 * Created on 3 giugno 2015, 23.20
 */

#ifndef CIFSAVER_H
#define	CIFSAVER_H

// Includes:
#include <string>
#include <iostream>

#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <LigandSet.h>
#include <Ligand.h>
#include <Saver.h>
#include <Debug.h>
#include <Protein.h>

#include "CifStructure.h"

// Global constants, typedefs, etc. (to avoid):

namespace Victor {
    namespace Biopool {

        /**
         * @brief Saves components (Atoms, Groups, etc.) in standard PDB format
         * */
        class CifSaver : public Saver {
        public:
            // CONSTRUCTORS/DESTRUCTOR:

            /**
             * Basic constructor. By default it writes sequence,
             * secondary structure and the term line. 
             * @param _output the output file object
             */
            CifSaver(ostream& _output = cout);

            // this class uses the implicit copy operator.

            virtual ~CifSaver();

            // PREDICATES:

            void endFile();

            // MODIFIERS:

            void setWriteSecondaryStructure();
            void setDoNotWriteSecondaryStructure();
            void setWriteSeqRes();
            void setDoNotWriteSeqRes();
            void setWriteAtomOnly();
            void setWriteAll();
            void setChain(char _ch);

            /**
             * Saves a group in CIF format.
             * @param group reference 
             * @return void
             */
            virtual void saveGroup(Group& gr);
            virtual void saveSideChain(SideChain& sc);
            virtual void saveAminoAcid(AminoAcid& aa);
            virtual void saveSpacer(Spacer& sp);
            virtual void saveLigand(Ligand& l);
            virtual void saveLigandSet(LigandSet& l);
            virtual void saveProtein(Protein& prot);

        protected:

        private:
            // HELPERS:
            
            // writes SEQRES entry
            void writeSeqRes(Spacer& sp);
            
            // writes secondary entries (SHEET, HELIX, etc.)
            void writeSecondary(Spacer& sp);
            
            // ATTRIBUTES 
            ostream& output; // output stream
            bool writeSeq, writeSecStr, writeTer;
	    
            // offsets that determine at which atom,
            // aminoacid and ligand number to start
            unsigned int atomOffset, ligandOffset;
            int aminoOffset;
            char chain; // chain ID
	    
	    bool atomGroupPrinted;
	    
	    CifStructure* cif;
        };

        inline void CifSaver::endFile() {
            output << "END\n";
        }

        inline void CifSaver::setWriteSecondaryStructure() {
            writeSecStr = true;
        }

        inline void CifSaver::setDoNotWriteSecondaryStructure() {
            writeSecStr = false;
        }

        inline void CifSaver::setWriteSeqRes() {
            writeSeq = true;
        }

        inline  void CifSaver::setDoNotWriteSeqRes() {
            writeSeq = false;
        }

        inline void CifSaver::setWriteAtomOnly() {
            writeSecStr = false;
            writeSeq = false;
            writeTer = false;
        }

        inline void CifSaver::setWriteAll() {
            writeSecStr = true;
            writeSeq = true;
            writeTer = true;
        }

        inline void CifSaver::setChain(char _ch) {
            chain = _ch;
        }

    }
} //namespace

#endif	/* CIFSAVER_H */


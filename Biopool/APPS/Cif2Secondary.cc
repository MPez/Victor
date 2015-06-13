/* 
 * File:   cif2secondary.cc
 * Author: marco
 *
 * Created on 12 giugno 2015, 17.42
 */

#include <string>
#include <GetArg.h>
#include <Spacer.h>
#include <CifLoader.h>
#include <IoTools.h>

using namespace Victor;
using namespace Victor::Biopool;

void sShowHelp() {
    cout << "Pdb 2 Secondary Structure converter\n"
	    << "\t H = helix, \t E = extended (strand, sheet), \t . = other.\n"
	    << "   Options: \n"
	    << "\t-i <filename> \t\t Input file for PDB structure\n"
	    << "\n";
}

int main(int argc, char** argv) {

    if (getArg("h", argc, argv)) {
	sShowHelp();
	return 1;
    };
    vector<char> allCh;
    string chainID = "!";
    string inputFile;
    getArg("i", inputFile, argc, argv, "!");

    if (inputFile == "!") {
	cout << "Missing file specification. Aborting. (-h for help)" << endl;
	return -1;
    }
    
    ifstream inFile(inputFile.c_str());
    if (!inFile) {
	ERROR("File not found.", exception);
    }

    CifLoader il(inFile);
    il.setNoHAtoms();
    allCh = il.getAllChains();
    
    for (unsigned int i = 0; i < allCh.size(); i++) {
	cout << "\t," << allCh[i] << ",";
    }
    cout << "\n";

    /*check on validity of chain: 
    if user select a chain then check validity
     else select first valid one by default*/
    if (chainID != "!") {
	bool validChain = false;
	for (unsigned int i = 0; i < allCh.size(); i++) {
	    if (allCh[i] == chainID[0]) {
		il.setChain(chainID[0]);
		cout << "Loading chain " << chainID << "\n";
		validChain = true;
		break;
	    }
	}
	if (!validChain) {
	    cout << "Chain " << chainID << " is not available\n";
	    return -1;
	}

    } else {
	chainID[0] = allCh[0];
	cout << "Using chain " << chainID << "\n";
    }
    
    Protein prot;
    prot.load(il);
    Spacer *sp;
    sp = prot.getSpacer(chainID[0]);

    allCh = il.getAllChains();
    cout << ">" << inputFile << "\n";
    
    for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
	switch (sp->getAmino(i).getState()) {
	    case HELIX:
		cout << "H";
		break;
	    case STRAND:
		cout << "E";
		break;
	    default:
		cout << ".";
	};
	if ((i + 1) % 60 == 0)
	    cout << "\n";
    }
    cout << "\n";

    return 0;
}


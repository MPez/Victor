/* 
 * File:   CifEditor.cc
 * Author: marco
 *
 * Created on 8 giugno 2015, 12.32
 */

#include <CifLoader.h>
#include <CifSaver.h>

using namespace Victor;
using namespace Victor::Biopool;

int main(int argc, char** argv) {
    if (argc != 3) {
        cout << "CIF Editor $Revision: 0.1 $ -- allows sequential manipulation of "
                << "protein structure backbone torsion angles" << endl;
        cout << "  Usage: \t\t CifEditor <input_filename> <output_filename> \n";
        return 1;
    }
    
    ifstream inFile(argv[1]);
    if (!inFile)
        ERROR("Input file not found.", exception);

    CifLoader cl(inFile);
    Protein prot;
    prot.load(cl);
    unsigned int zero = 0;
    Spacer sp = *(prot.getSpacer(zero));

    cout << "Editing " << argv[1] << " output goes to " << argv[2] << "\n";

    int aaid = -1;
    do {
        cout << "Aminoacid# or -1: ";
        cin >> aaid;
        if (aaid <= -1) {
            cout << "Bye.\n";
            return 0;
        }

        if (aaid >= (int) sp.sizeAmino()) {
            cout << "\t Invalid aa#!\n";
        } else {
            double newVal = 999;
            cout << "\t " << aaid << "  " << sp.getAmino(aaid).getType() << "\n";
            cout << "Phi=   " << sp.getAmino(aaid).getPhi()
                    << "\t new phi or 999:   ";
            cin >> newVal;
            if (newVal != 999)
                sp.getAmino(aaid).setPhi(newVal);
            cout << "Psi=   " << sp.getAmino(aaid).getPsi()
                    << "\t new psi or 999:   ";
            cin >> newVal;
            if (newVal != 999)
                sp.getAmino(aaid).setPsi(newVal);
            cout << "Omega= " << sp.getAmino(aaid).getOmega()
                    << "\t new omega or 999: ";
            cin >> newVal;
            if (newVal != 999)
                sp.getAmino(aaid).setOmega(newVal);

            ofstream outFile2(argv[2]);

            if (!outFile2)
                ERROR("Couldn't write output file.", exception);

            CifSaver pss2(outFile2);

            sp.save(pss2);
        }
    } while (aaid != -1);

    return 0;
}


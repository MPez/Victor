/* 
 * File:   CifCorrector.cc
 * Author: marco
 *
 * Created on 12 giugno 2015, 18.48
 */

#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <CifLoader.h>
#include <CifSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Victor;
using namespace Victor::Biopool;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc != 2) {
        cout << "Cif Corrector $Revision: 0.1 $ -- adds missing oxygen atoms to "
                << "protein structure backbones" << endl;
        cout << "  Usage: \t\t CifCorrector <filename> \n";
        return 1;
    };

    Spacer sp;
    ifstream inFile(argv[1]);
    if (!inFile)
        ERROR("File not found.", exception);

    CifLoader il(inFile);
    sp.load(il);

    for (unsigned int i = 0; i < sp.sizeAmino(); i++)
        sp.getAmino(i).addMissingO();

    ofstream outFile2(argv[1]);

    if (!outFile2)
        ERROR("Couldn't write file.", exception);

    CifSaver pss2(outFile2);
    sp.save(pss2);
    
    return 0;
}


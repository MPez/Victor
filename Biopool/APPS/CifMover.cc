/* 
 * File:   cifMover.cc
 * Author: marco
 *
 * Created on 12 giugno 2015, 18.39
 */

#include <string>
#include <GetArg.h>
#include <CifLoader.h>
#include <CifSaver.h>
#include <vector3.h>
#include <matrix3.h>
#include <IntCoordConverter.h>

using namespace Victor;
using namespace Victor::Biopool;

// minimum distance between neighbouring CAs
const double LAMBDA = 1.5;

void sShowHelp() {
    cout << "CIF Shifter\n"
	    << "Allows to move all residues in file by fixed offset.\n"
	    << " Options: \n"
	    << "\t-i <filename> \t\t Input CIF file\n"
	    << "\t-o <filename> \t\t Output CIF file\n"
	    << "\t[-r <repeat>] \t\t Repeat unit length (default = no rotation)\n"
	    << "\t[-l <double>] \t\t Angle lambda factor (default = 0.1)\n"
	    << "\t[-s <number>] \t\t Start residue of fragment (default = first)\n"
	    << "\t[-e <number>] \t\t End residue of fragment (default = last)\n"
	    << "\n";
}

void sAddLine() {
    cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

int main(int argc, char** argv) {

    // --------------------------------------------------
    // 0. treat options
    // --------------------------------------------------

    if (getArg("h", argc, argv)) {
	sShowHelp();
	return 1;
    };

    string inputFile, outputFile;
    unsigned int startOffset, endOffset, repeatLength;
    double lambdaAngle;

    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");
    getArg("s", startOffset, argc, argv, 0);
    getArg("e", endOffset, argc, argv, 9999);
    getArg("r", repeatLength, argc, argv, 9999);
    getArg("l", lambdaAngle, argc, argv, 0.1);

    vgVector3<double> transOff;
    for (unsigned int i = 0; i < 3; i++)
	transOff[i] = 0.0;

    if ((inputFile == "!") || (outputFile == "!")) {
	cout << "Missing file specification. Aborting. (-h for help)" << endl;
	return -1;
    }

    // --------------------------------------------------
    // 1. read structure
    // --------------------------------------------------

    ifstream inFile(inputFile.c_str());

    if (!inFile)
	ERROR("File does not exist.\n", exception);
    
    CifLoader cl(inFile);
    Protein prot;
    prot.load(cl);
    unsigned int zero = 0;
    Spacer sp = *(prot.getSpacer(zero));
    
    inFile.close();


    endOffset = sp.getIndexFromPdbNumber(endOffset);
    if (startOffset > 0)
	startOffset = sp.getIndexFromPdbNumber(startOffset);

    // --------------------------------------------------
    // 2. rotate spacer
    // --------------------------------------------------

    if (repeatLength < 9999) {
	IntCoordConverter icc;

	// find offset
	vgVector3<double> firstA = sp.getAmino(endOffset
		- (3 * repeatLength / 4))[CA].getCoords()
		- sp.getAmino(endOffset)[CA].getCoords();

	vgVector3<double> firstB = sp.getAmino(endOffset
		- (1 * repeatLength / 4))[CA].getCoords()
		- sp.getAmino(endOffset)[CA].getCoords();

	vgVector3<double> firstNorm = (firstA.normalize()).cross(firstB.normalize());

	vgVector3<double> secondA = sp.getAmino(endOffset + repeatLength
		- (3 * repeatLength / 4))[CA].getCoords()
		- sp.getAmino(endOffset)[CA].getCoords();

	vgVector3<double> secondB = sp.getAmino(endOffset + repeatLength
		- (1 * repeatLength / 4))[CA].getCoords()
		- sp.getAmino(endOffset)[CA].getCoords();

	vgVector3<double> secondNorm = (secondA.normalize()).cross(secondB.normalize());

	firstNorm.normalize();
	secondNorm.normalize();

	double scalar = icc.getAngle(firstNorm, secondNorm) * lambdaAngle;

	cout << "Scalar = " << setw(5) << setprecision(3)
		<< RAD2DEG * scalar / lambdaAngle << "\n";


	vgVector3<double> axis = (firstNorm.normalize()).cross(secondNorm.normalize());
	vgMatrix3<double> res = vgMatrix3<double>::createRotationMatrix(axis, scalar);

	sp.getAmino(startOffset)[N].addRot(res);
    }

    // --------------------------------------------------
    // 3. translate spacer
    // --------------------------------------------------

    if (endOffset < sp.sizeAmino() - 1) {
	sp.getAmino(endOffset).unbindOut(sp.getAmino(endOffset + 1));
	sp.getAmino(endOffset)[C].unbindOut(sp.getAmino(endOffset + 1)[N]);

	// find offset
	vgVector3<double> first = sp.getAmino(endOffset)[C].getCoords();
	vgVector3<double> second = sp.getAmino(endOffset + 1)[N].getCoords();

	double d = sp.getAmino(endOffset)[C].distance(
		sp.getAmino(endOffset + 1)[N]);

	double frac = (d - LAMBDA) / d;

	for (unsigned int i = 0; i < 3; i++)
	    transOff[i] = (second[i] - first[i]) * frac;
    } else {
	if (startOffset != 0)
	    ERROR("Both start and end offset are undefined.", exception);

	endOffset = sp.sizeAmino() - 1;

	// NB: unbind has to be reversed after moving the atoms if the model
	// is to be used further in the same program
	sp.getAmino(startOffset - 1).unbindOut(sp.getAmino(startOffset));
	sp.getAmino(startOffset - 1)[C].unbindOut(sp.getAmino(startOffset)[N]);

	// find offset
	vgVector3<double> first = sp.getAmino(startOffset + 1)[N].getCoords();
	vgVector3<double> second = sp.getAmino(startOffset)[C].getCoords();

	double d = sp.getAmino(startOffset + 1)[N].distance(
		sp.getAmino(startOffset)[C]);

	double frac = (d - LAMBDA) / d;

	for (unsigned int i = 0; i < 3; i++)
	    transOff[i] = (second[i] - first[i]) * frac;
    }

    sp.getAmino(startOffset)[N].addTrans(transOff);

    // --------------------------------------------------
    // 4. write model to disk
    // --------------------------------------------------

    ofstream outFile(outputFile.c_str());
    if (!outFile) {
	ERROR("File not found.", exception);
    }
    CifSaver cs(outFile);
    sp.save(cs);
    outFile.close();

    return 0;
}


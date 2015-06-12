/* 
 * File:   TestCifStructure.h
 * Author: marco
 *
 * Created on 10 giugno 2015, 10.37
 */

#ifndef TESTCIFSTRUCTURE_H
#define	TESTCIFSTRUCTURE_H

#include <iostream>
#include <fstream>

#include <cppunit/TestFixture.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>

#include <IoTools.h>
#include <String2Number.h>
#include <CifStructure.h>

using namespace std;
using namespace Victor::Biopool;
using namespace CppUnit;

class TestCifStructure : public TestFixture {
public:

    TestCifStructure() {
    }

    virtual ~TestCifStructure() {
    }

    static Test* suite() {
	TestSuite* suiteOfTests = new TestSuite("TestCifLoader");

	suiteOfTests->addTest(new TestCaller<TestCifStructure>("Get group column number",
		&TestCifStructure::testGetGroupColumnNumber));

	suiteOfTests->addTest(new TestCaller<TestCifStructure>("Get group field",
		&TestCifStructure::testGetGroupField));
	
	suiteOfTests->addTest(new TestCaller<TestCifStructure>("Get inline field",
		&TestCifStructure::testGetInlineField));

	return suiteOfTests;
    }

    void setUp() {
	// inizialize CifLoader
	string path = getenv("VICTOR_ROOT");
	string input = path + "Biopool/Tests/data/modelTest.cif";
	inFile = new ifstream(input.c_str());
	cs = new CifStructure(*inFile);
	
	// initialize test parameters
	idColumn = 1;
	atomId = 1;
	headerLine = "_entry.id   25C8 ";
	headerId = "25C8";
    }

    void tearDown() {
	delete cs;
	delete inFile;
    }

private:

    void testGetGroupField() {
	string line = readLine(cs->getInput());
	cs->parseGroup("atom", line);
	int field = stoiDEF(cs->getGroupField("atom", line, 
		cs->getGroupColumnNumber("atom", "atom id")));
	
	CPPUNIT_ASSERT_EQUAL(atomId, field);
    }

    void testGetGroupColumnNumber() {
	string line = readLine(cs->getInput());
	cs->parseGroup("atom", line);
	int col = cs->getGroupColumnNumber("atom", "atom id");
	
	CPPUNIT_ASSERT_EQUAL(idColumn, col);
    }
    
    void testGetInlineField() {
	string field = cs->getInlineField(headerLine);
	
	CPPUNIT_ASSERT_EQUAL(headerId, field);
    }
    
    int idColumn;
    int atomId;
    string headerLine;
    string headerId;
    
    ifstream* inFile;
    CifStructure* cs;
    
};

#endif	/* TESTCIFSTRUCTURE_H */


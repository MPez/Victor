/*
 * File:   TestCifLoader.h
 * Author: marco
 *
 * Created on 9-giu-2015, 10.45.39
 */

#ifndef TESTCIFLOADER_H
#define	TESTCIFLOADER_H

#include <iostream>
#include <fstream>

#include <cppunit/TestFixture.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>

#include <CifLoader.h>

using namespace std;
using namespace Victor::Biopool;
using namespace CppUnit;

class TestCifLoader : public TestFixture {
public:

    TestCifLoader() {
    }

    virtual ~TestCifLoader() {
    }

    static Test* suite() {
	TestSuite* suiteOfTests = new TestSuite("TestCifLoader");

	suiteOfTests->addTest(new TestCaller<TestCifLoader>("Test 1 - Get max numebers of models",
		&TestCifLoader::testGetMaxModels));

	suiteOfTests->addTest(new TestCaller<TestCifLoader>("Test 2 - Get all chain ids",
		&TestCifLoader::testGetAllChains));

	return suiteOfTests;
    }

    void setUp() {
	// inizialize CifLoader
	string path = getenv("VICTOR_ROOT");
	string input = path + "Biopool/Tests/data/modelTest.cif";
	inFile = new ifstream(input.c_str());
	cl = new CifLoader(*inFile);
	
	// initialize test parameters
	maxModel = 5;
	chainIds.push_back('A');
	chainIds.push_back('B');
	chainIds.push_back('C');
    }

    void tearDown() {
	delete cl;
	delete inFile;
    }

private:

    void testGetMaxModels() {
	int max = cl->getMaxModels();
	CPPUNIT_ASSERT_EQUAL(max, maxModel);
    }

    void testGetAllChains() {
	vector<char> testChain = cl->getAllChains();
	CPPUNIT_ASSERT(chainIds[0] == testChain[0] &&
		chainIds[1] == testChain[1] &&
		chainIds[2] == testChain[2]);
    }

    int maxModel;
    vector<char> chainIds;
    
    ifstream* inFile;
    CifLoader* cl;
};


/*
#include <cppunit/extensions/HelperMacros.h>
#include "CifLoader.h"
#include <fstream>

using namespace::Victor::Biopool;
using namespace::CppUnit;

class TestCifLoader : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestCifLoader);

    CPPUNIT_TEST(testGetMaxModels);
    CPPUNIT_TEST(testGetAllChains);

    CPPUNIT_TEST_SUITE_END();

public:
    TestCifLoader();
    virtual ~TestCifLoader();
    void setUp();
    void tearDown();

private:
    void testGetMaxModels();
    void testGetAllChains();
    
    CifLoader* cl;
    int maxModel;
    vector<char> chainIds;
};
 */
#endif	/* TESTCIFLOADER_H */


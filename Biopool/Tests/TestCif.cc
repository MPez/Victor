/*
 * File:   TestCif.cc
 * Author: marco
 *
 * Created on 9-giu-2015, 10.45.40
 */

/*
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
*/

#include <iostream>
#include <cppunit/ui/text/TestRunner.h>


#include <TestCifLoader.h>

using std::cout;
using std::endl;

int main() {
    
    CppUnit::TextUi::TestRunner runner;

    cout << "Creating Test Suites:" << endl;
    
    runner.addTest(TestCifLoader::suite());
    
    cout << "Running the unit tests." << endl;
    
    runner.run();

    return 0;
    
    /*
    // Create the event manager and test controller
    CPPUNIT_NS::TestResult controller;

    // Add a listener that colllects test result
    CPPUNIT_NS::TestResultCollector result;
    controller.addListener(&result);

    // Add a listener that print dots as test run.
    CPPUNIT_NS::BriefTestProgressListener progress;
    controller.addListener(&progress);

    // Add the top suite to the test runner
    CPPUNIT_NS::TestRunner runner;
    runner.addTest(CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest());
    runner.run(controller);

    // Print test in a compiler compatible format.
    CPPUNIT_NS::CompilerOutputter outputter(&result, CPPUNIT_NS::stdCOut());
    outputter.write();

    return result.wasSuccessful() ? 0 : 1;
    */
}

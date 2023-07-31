Getting started with Software Testing 
=====================================


(To do: It should be clear from the title of the page what is the reading order of this page in the developer documentation. There is a quickstart page as well. To begin with, we can specify the order in index.rst)

Alternate script: An effective way to walk through a newbie to software tests is to let them dive straight into the practical aspects of how they are implemented and and then explain the theory behind it. This is the bottom up approach taken in this page. Top down approach can disconnect the reader as they will see only theory and not connect it with the tests that are implemented in the code. Write this page as an interactive tutorial which tells a story and engages the reader.  

Opening paragraph: This page is guide for the developers of aeolis on software tests. It explains what software tests are, why they are important and how they are implemented in aeolis. We take the bottom up approach to show you what software tests really are and in which form they exist in aeolis. This is meant to give a more practical feel of how tests look like in practice. 

- Tests are implemented using the pytest (add URL) framework. To run the tests run 'pytest' from the root of the repository.

You will see the below prompt on the command line. 
"test session starts..."

(Add screenshot.)

What is happening here is that by executing pytest on the command line you have started the execution of a set of scripts that are written to test the software. These scripts are called tests. You can find these scripts in aeolis/tests/ folder.

**So, software tests are actually code that is written to test the software. Software tests are essentially a piece of code, for example, Python scripts, that are run to check whether a target software, for example AeoLiS, is working as expected or not. The tests generate a Pass or Fail status as the output upon execution to indicate this.** 


- If you want to see what the scripts are executing, you can abort the current operation with Ctrl+C and rerun pytest with the -s option. 

- After all the scripts are run, you see the below output on the command line. 

```=========================== 2 passed in 0.04s ===========================```

This indicates that all the tests have passed. If tests fail, it will be indicated in the output.


Contents of a software test/What's inside a software test?
----------------------------------------------------------

At this point, you can peek into one of the test scripts and see what it looks like. For example, cat aeolis/tests/functional_tests/test_aeolis.py. **The test script is similar to a normal Python module with a collection of functions, classes, etc.** The only difference is that the functions in the test script are prefixed with the word `test_`. This enables them to be discovered and run by pytest upon the execution of the command pytests on the command line. 


Writing tests
-------------

In a nutshell, the test scripts are a collection of functions that are prefixed with the word `test_`. These functions call the functions in various modules within the aeolis source code with certain inputs and check whether the output is as expected. If the output is as expected, the test passes. If the output is not as expected, the test fails. This is the basic idea behind software tests. For an example, see the section Example formatting log messages in unit testing page to see how the test functions are written.


Automated testing 
-----------------
We ran tests manually now, but they can be configured to run automatically when the code is changed. This is called continuous integration. This is done using the GitHub Actions framework. The configuration for this is in the file .github/workflows/python-package.yml. More details are in the automated testing page. **Thus, by running software tests alongside the code development process, we can ensure that the code is always working as expected.** 



Types of tests
--------------

Software tests can be classified into different types based on the scope of the software that they test:

- Unit tests: test a single function or a small piece of code. See the page unit testing to see the unit tests in AeoLiS.

- Integration tests: test how different pieces of code work together    

- Functional tests: Functional testing is a type of black-box testing where the software as a whole is tested by feeding it input and examining the output, and internal program structure is rarely considered. Functional testing usually describes what the system does. Functional testing does not imply that you are testing a function (method) of a module or class. Functional testing tests a slice of functionality of the whole system. 

- End-to-end tests: test the entire software


Software tests in aeolis
------------------------

aeolis currently has two types of tests:

- Unit tests: test the functions in the `aeolis` package. These are located in aeolis/tests folder.
- Functional tests: treats the software as a black box and tests the functionality of the software as a whole by checking its output against a particular input. These are located in aeolis/tests/functional_tests folder.

Read the implementation details of the above tests in aeolis on the pages [Unit testing](unit_testing.md) and [Functional testing](functional_testing.md).


=====================================
Getting started with Software Testing 
=====================================

This document aims to introduce aeolis developers to the fundamental concepts of software testing. By the end of this document, you will be able to answer the following questions: 

- What are software tests?
- What is inside a software test? What do software tests look like?
- How do software tests work?
- What do software tests do?

To answer these questions, we will use the existing tests in aeolis as the starting point. 

You will need the aeolis source code and the aeolis development environment set up on your machine to follow this document. If you do not have the development environment set up, please follow the instructions in the [Development environment setup](development_environment_setup.md) page.


Software tests in aeolis
------------------------

In order to know what do software tests look like and how they work in practice, we start by running the existing tests in aeolis and observing the results.

Run 'pytest' from the root of the repsoitory to run the existing tests in aeolis. 

```pytest```

Pytest is a Python package that is used to run software tests. It is a framework that is used to write and run software tests. It is a third-party package that is not part of the Python standard library. It is installed as a dependency when you install aeolis. 


You will see the below prompt on the command line. 
"test session starts..."

(Add screenshot.)


After all the scripts are run (~ 4 minutes), you see the below output on the command line. 

```=========================== 2 passed in 0.04s ===========================```

Add screenshot

This output indicates that all the tests have passed. If tests fail, it will be indicated in the output.


Key observations
^^^^^^^^^^^^^^^^

- Executing `pytest` on the command line starts the execution of a set of Python scripts. These scripts start with the name `test_`. These scripts are called tests. You can find these scripts in aeolis/tests/ folder.

- **Software tests are essentially a piece of code, for example, Python scripts, that are run to check whether a target software, for example Aeolis, is working as expected or not.** 

- **Tests produce a pass or fail status as the output which indicates whether the target software is working as expected or not.**



What's inside a software test? What does it look like?
------------------------------------------------------

At this point, you can peek into one of the test scripts and see what it looks like. For example, cat aeolis/tests/functional_tests/test_aeolis.py. **The test script is similar to a normal Python module with a collection of functions, classes, etc.** The only difference is that the functions in the test script are prefixed with the word `test_`. This enables them to be discovered and run by pytest upon the execution of the command pytests on the command line. 


If you want to see what the scripts are executing, you can abort the current operation with Ctrl+C and rerun pytest with the -s option. This will print the output of the scripts to the command line. 

```pytest -s```


Key takeaways
^^^^^^^^^^^^^

- In a nutshell, a test script is a collection of functions that are prefixed with the word `test_`. These functions call the functions in various modules within the aeolis source code with certain inputs and check whether the output is as expected. If the output is as expected, the test passes. If the output is not as expected, the test fails. This is the basic idea behind software tests. For an example, see the section Example formatting log messages in unit testing page to see how the test functions are written.


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


aeolis currently has two types of tests: unit tests and functional tests.

Read the implementation details of the above tests in aeolis in the document [Unit testing](unit_testing.md) and [Functional testing](functional_testing.md).


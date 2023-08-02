=====================================
Getting started with Software Testing 
=====================================

This document aims to introduce aeolis developers to the fundamental concepts of software testing. By the end of this document, you will be able to answer the following questions: 

- What are software tests?
- What do software tests do?
- How do software tests work?

Existing tests in aeolis are used as the starting point to answer the above questions.

You will need the aeolis source code and the aeolis development environment set up on your machine to follow this document. If you do not have the development environment set up, please follow the  instructions in the section :doc:`Set up <quickstart>`.

*In order to know what software tests are and how they work in practice, follow the below steps to run the existing tests in aeolis and observe the results.*

1. Clone the aeolis-python GitHub repository on your local machine and execute ```pytest``` from the root of the repsoitory. This produces the below output on the console. 

.. code-block:: console

    :~/aeolis-python$ pytest
    ============================================================================ test session starts ============================================================================
    platform linux -- Python 3.8.10, pytest-7.2.2, pluggy-1.0.0
    rootdir: /~/aeolis-python
    plugins: cov-4.0.0
    collected 18 items

    aeolis/tests/test_model.py ........                                                                                                                                   [ 44%]
    aeolis/tests/test_utils.py ........                                                                                                                                   [ 88%]
    aeolis/tests/regression_tests/test_netCDF_file_content.py


- Executing `pytest` on the command line starts the execution of a set of Python scripts that start with the word ```test_``` and are located in the aeolis/tests/ folder. These scripts are the software tests that are used to test the aeolis source code.

- Pytest is a Python package and testing framework that is used to run software tests. It is a third-party package that is not part of the Python standard library. It is installed as a dependency when you install aeolis. It provides a command line interface to run the tests. Executing `pytest` on the command line runs all the tests in the aeolis/tests/ folder.

- It currently takes approximately 2.5 minutes for all the test scripts to run. Once the tests finish running, you will see the below output on the console.



.. code-block:: console

    :~/aeolis-python$ pytest
    ============================================================================ test session starts ============================================================================
    platform linux -- Python 3.8.10, pytest-7.2.2, pluggy-1.0.0
    rootdir: /~/aeolis-python
    plugins: cov-4.0.0
    collected 18 items

    aeolis/tests/test_model.py ........                                                                                                                                   [ 44%]
    aeolis/tests/test_utils.py ........                                                                                                                                   [ 88%]
    aeolis/tests/regression_tests/test_netCDF_file_content.py .                                                                                                           [ 94%]
    aeolis/tests/regression_tests/test_netCDF_file_creation.py .                                                                                                          [100%]

    ============================================================================= warnings summary ==============================================================================
    aeolis/tests/regression_tests/test_netCDF_file_content.py: 4527 warnings
    aeolis/tests/regression_tests/test_netCDF_file_creation.py: 4527 warnings
    /~/venvs/aeolis-dev/lib/python3.8/site-packages/numpy/matrixlib/defmatrix.py:69: PendingDeprecationWarning: the matrix subclass is not the recommended way to 
        represent matrices or deal with linear algebra (see https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html). 
        Please adjust your code to use regular ndarray.
        return matrix(data, dtype=dtype, copy=False)

    -- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
    =============================================================== 18 passed, 9054 warnings in 118.43s (0:01:58) ===============================================================


- The test session section displays the status of the tests. It shows the number of tests that passed and the number of tests that failed. In the above example, all the tests passed.


2. To know what a software test exactly is and what it contains, open one of the test scripts in aeolis/tests/ folder and browse through its contents. For example, open aeolis/tests/regression_tests/test_aeolis.py.

   The test script is similar to a normal Python module with a collection of functions, classes, etc. The only difference is that the functions in the test script are prefixed with the word `test_`. This enables them to be discovered and run by pytest upon the execution of the command pytests on the command line. 


Software tests in a nutshell
----------------------------

Based on the observations after running the tests in aeolis and examining the test scripts, below are the key points to remember about software tests:

- A software test is essentially a piece of code, such as Python scripts, that is executed to verify whether a target software, for example, aeolis, is functioning as expected or not.

- Tests produce a pass or fail status as the output to indicate whether the target software is working as expected or not.

- A test script is a collection of functions that are prefixed with the word `test_`. These functions call the functions in various modules within the aeolis source code with certain inputs and check whether the output is as expected. If the output is as expected, the test passes. If the output is not as expected, the test fails. This is the basic idea behind software tests. For an example, see the section *Example: formatting log messages* in :doc:`unit testing <unit-testing>` to learn how to write test functions.


Types of tests
--------------

Software tests can be classified into different types based on the scope of the software that they test:

- Unit tests: Test a single function or a small piece of code.

- Integration tests: Test how different pieces of code work together.     

- Regression tests: Regression testing is a type of black-box testing where the software as a whole is tested by feeding it input and examining the output, and internal program structure is rarely considered. regression testing usually describes what the system does. regression testing does not imply that you are testing a function (method) of a module or class. regression testing tests a slice of regressionity of the whole system. 


Software testing setup in aeolis
--------------------------------

- aeolis currently has two types of tests: unit tests and regression tests. Read the implementation details of the above tests in aeolis in the document :doc:`unit testing <unit-testing>` and :doc:`regression testing <regression-testing>`.
- aeolis-python GitHub repository has automated testing enabled to run tests automatically on every push to main branch to catch any potential bugs that the new code may introduce. Read more about this in the document :doc:`automated testing <automated-testing>`.




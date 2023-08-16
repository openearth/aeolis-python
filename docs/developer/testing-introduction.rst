=====================================
Getting started with Software Testing 
=====================================

This document introduces AeoLiS developers to the fundamental concepts of software testing using `Pytest <https://docs.pytest.org/en/7.1.x/getting-started.html>`_. The following topics are included: 

- What are software tests?
- What do software tests do?
- How software tests work?

We show some examples using existing tests in AeoLiS. You will need the AeoLiS source code and the development environment if you want to reproduce the examples. Follow the instruction in the :doc:`quick start <quickstart>` to set up a development environment.


What are software tests?
------------------------

To show what software tests are and how they work, you will follow some examples using Pytest.

1. Clone the aeolis-python GitHub repository to your machine and execute `pytest` from the root of the repsoitory using the terminal. This produces following output: 

.. code-block:: console

    :~/aeolis-python$ pytest
    ==================================== test session starts ==================================== 
    platform linux -- Python 3.8.10, pytest-7.2.2, pluggy-1.0.0
    rootdir: /~/aeolis-python
    plugins: cov-4.0.0
    collected 18 items

    aeolis/tests/test_model.py ........                                                    [ 44%]
    aeolis/tests/test_utils.py ........                                                    [ 88%]
    aeolis/tests/regression_tests/test_netCDF_file_content.py

.. note:: 

   It takes approximately 2.5 minutes for all the test files to run. Once the tests finish running


- Executing `pytest` on the command line starts the execution of a set of Python files which names start with *test_* and are located in the `aeolis/tests/`` folder. These files contain code that it is used to test the aeolis source code.

- Pytest is a Python package and testing framework that is used to test software. It is a third-party package that it is installed when you install aeolis using the development depencies. 

- Once the tests finish running, you will see the below output on the console.


.. code-block:: console

    :~/aeolis-python$ pytest
    ==================================== test session starts ===================================
    platform linux -- Python 3.8.10, pytest-7.2.2, pluggy-1.0.0
    rootdir: /~/aeolis-python
    plugins: cov-4.0.0
    collected 18 items

    aeolis/tests/test_model.py ........                                                   [ 44%]
    aeolis/tests/test_utils.py ........                                                   [ 88%]
    aeolis/tests/regression_tests/test_netCDF_file_content.py .                           [ 94%]
    aeolis/tests/regression_tests/test_netCDF_file_creation.py .                          [100%]

    ====================================== warnings summary ==================================== 
    aeolis/tests/regression_tests/test_netCDF_file_content.py: 4527 warnings
    aeolis/tests/regression_tests/test_netCDF_file_creation.py: 4527 warnings
    /~/venvs/aeolis-dev/lib/python3.8/site-packages/numpy/matrixlib/defmatrix.py:69: 
        PendingDeprecationWarning: the matrix subclass is not the recommended way to 
        represent matrices or deal with linear algebra 
        (see https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html). 
        Please adjust your code to use regular ndarray.
        return matrix(data, dtype=dtype, copy=False)

    -- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
   ======================= 18 passed, 9054 warnings in 118.43s (0:01:58) ======================


- The test **session section** displays the status of the tests. It shows the number of tests that passed and the number of tests that failed. In this example, all the tests passed. The **warnings summary** section displays the number of warnings that were generated during the execution of the tests. Warnings are a feature of Pytest that checks for any potential issues in the code, but do not affect the result of the tests.

[Continue here]

1. To know what a software test exactly is and what it contains, open one of the test files in aeolis/tests/ folder and browse through its contents. For example, open aeolis/tests/regression_tests/test_aeolis.py.

   The test script is similar to a normal Python module with a collection of functions, classes, etc. The only difference is that the functions in the test script are prefixed with the word `test_`. This enables them to be discovered and run by pytest upon the execution of the command pytests on the command line. 


Software tests in a nutshell
----------------------------

Based on the observations after running the tests in aeolis and examining the test files, below are the key points to remember about software tests:

- A software test is essentially a piece of code, such as Python files, that is executed to verify whether a target software, for example, aeolis, is functioning as expected or not.

- Tests produce a pass or fail status as the output to indicate whether the target software is working as expected or not.

- A test script is a collection of functions that are prefixed with the word `test_`. These functions call the functions in various modules within the aeolis source code with certain inputs and check whether the output is as expected. If the output is as expected, the test passes. If the output is not as expected, the test fails. This is the basic idea behind software tests. For an example, see the section *Example: formatting log messages* in :doc:`unit testing <unit-testing>` to learn how to write test functions.


Types of tests
--------------

Software tests can be classified into different types based on the scope of the software that they test:

- Unit tests: Test a single function or a small piece of code. To learn how to write unit tests for aeolis, read the section :doc:`unit testing <unit-testing>`.


- Regression tests: Regression testing is a type of black-box testing where the software as a whole is tested by feeding it input and examining the output, and internal program structure is rarely considered. 

  AeoLiS currently has the following two regression tests:

    - `test_netCDF_content.py <https://github.com/openearth/aeolis-python/blob/main/aeolis/tests/integration_tests/test_netCDF_file_content.py>`_: tests whether the content of the netCDF files generated by the latest version of AeoLiS for a given input file is the same as the content of the netCDF files generated by the previous version of AeoLiS for the same input file.

    - `test_netCDF_creation.py <https://github.com/openearth/aeolis-python/blob/main/aeolis/tests/integration_tests/test_netCDF_file_creation.py>`_: Tests whether a netCDF file is created upon a successful execution of the model.


- Integration tests: Test how different pieces of code work together.     


Automated testing in aeolis
---------------------------

- To catch any bugs introduced by new code changes, the test suite in Aeolis (unit tests + regression tests) runs automatically on every push to the main branch. This automated testing is enabled in the aeolis repository using `GitHub Actions <https://docs.github.com/en/actions/automating-builds-and-tests/building-and-testing-python>`_. If the tests fail, the new changes are not pushed to main. The test output status is displayed using a green tick (pass) or a red cross (fail) next to the commit message. 

- The file `.github/workflows/python-app.yml <https://github.com/openearth/aeolis-python/blob/main/.github/workflows/python-app.yml>`_ contains the instructions for the automated testing. The tests are run using the pytest framework. The tests are run on a virtual machine (Ubuntu 20.04) using the GitHub-hosted runner and on Python versions 3.8-3.11. This workflow is configured using GitHub Actions to run on every push to the main branch. 

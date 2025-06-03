=====================================
Introduction to Software Testing 
=====================================

This section is an introduction to the fundamental concepts of software testing using `Pytest <https://docs.pytest.org/en/7.1.x/getting-started.html>`_. It provides a quick start to software testing and shows some examples using existing tests in AeoLiS. 


What is software testing?
-------------------------

In short, **software testing** is the process of verifying and validating that a *software does what it is supposed to do*. The collection of *actions* that a software does is called the *behaviour* of the software, and in order to write test for a software, the developer needs to understand its expected behaviour. The expected behaviour of a software is defined by the user's *requirements*. Thefore, writing useful tests for a software is easier when the developer has a clear understanding of the user's requirements.

Consider a software with a single Python function which compares two strings:

.. code-block:: python

    def compare_strings(string1: srt, string2: srt) -> bool:
        """Compare two strings"""
        if string1 == string2:
            return True
        else:
            return False


The expected behaviour of this function is that it should return ``True`` if the two strings are equal and ``False`` otherwise.  Now that we know the expected behaviour of this code, we  can write tests to verify whether the function is working as expected or not. For example, we can test the following cases: firts, we can write a test that calls the function with two equal strings and checks whether the output is ``True``. Second, we can write a test that calls the function with two different strings and checks whether the output is ``False``.

.. code-block:: python

    # Case 1: equal strings must return True
    compare_strings("aeolis", "aeolis") 

    # Case 2: different strings must return False
    compare_strings("aeolis", "aeolis-python") 

    # Case 3: 
    compare_strings("aeolis", "Aeolis")


What about a third case? For example, when the two strings contain the same letters but one of them is in uppercase and the other is in lowercase. Should the function return ``True`` or ``False``? This is a question that the developer needs to answer based on the user's requirements. If the user's requirements state that the function should return ``True`` in that case, then the developer can write a test to check that as well. In short, *it is not possible to write meaningful tests for a software without knowing its expected behaviour.* A way to find out which tests need to be written for a software is to ask: **How the software should behave in this particular case?**  The answers to such quesiton lay on understanding the user's requirements.

.. note:: 

   - The example above is somehow simplistic. In a real case, the expected behaviour of a software is not always clear, and the developer needs to discuss with the users and other developers to understand which behaviour is expected.
   - Keep in mind that writing tests for a software is an iterative process. As the software evolves, the tests need to be updated to reflect the new behaviour of the software. Software tests provide feedback to improve the source code; writing tests can reveal cases that the software did not consider before but are imperative to handle. As a conseguence, the developer may decide to update the source code to handle such cases. 


How software tests work?
------------------------

To demonstrate how software tests work we will use Pytest, the test suite used in AeoLiS. Pytest is a Python package and testing framework that is used to test software. Pytest is installed  if the development dependencies are installed: 

1. Follow the instruction in the :doc:`quick start <quickstart>` to set up a development environment.

2. Clone the *aeolis-python* repository to your machine and execute ``pytest`` from the root of the repsoitory using the terminal. This produces following output: 

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

.. important:: 

   It takes approximately 2.5 minutes for all the test to run. 
   Executing ``pytest`` on the command line starts the execution of a set of Python files which names start with **test_** and are located in the ``aeolis/tests/`` folder. These files contain code that it is used to test the aeolis source code.

- Once the tests finish running, you will see an output similar to this in the console:

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


The test **session section** displays the status of the tests. It shows the number of tests that passed and the number of tests that failed. In this example, all the tests have passed. The **warnings summary** section displays the number of warnings that were generated during the execution of the tests. Warnings are a feature of Pytest that checks for any potential issues in the code, but do not affect the result of the tests.

What are software tests?
------------------------

Software tests are pieces of code that verify whether a target software is functioning as expected or not. The target software can be a single function, a module, or a whole application. The output of a test is either a pass or fail status. If the output is as expected, the test passes. If the output is not as expected, the test fails.


1. To see some examples of software tests, open one of the test files in ``aeolis/tests/`` directory and browse through its content. For example, open ``aeolis/tests/test_utils.py``. The test files are similar to a normal Python module with a collection of functions, classes, etc. The only difference is that by convention class names start with *Test* and test function are prefixed with the word ``test_``. The later, enables  Pytest to discovered them, and execute them upon the running the command ``pytest`` on the console. 

2. Focus your attention of the ``TestIsIterable``. This class tests if the ``is_iterable`` function in ``aeolis/utils.py`` works correctly. The ``is_iterable`` function checks whether its input is a Python iterable or not. The expected behaviour of this function is that it should return ``True`` if the object is iterable and ``False`` otherwise. The class ``TestIsIterable`` contains two test functions: ``test_is_iterable`` and ``test_is_not_iterable``. The first test function calls the ``is_iterable`` function with an input that we know before hand it an iterable object (a numpy array in this case) and checks whether the output is ``True``. The second test function calls the ``is_iterable`` function with a known-non-iterable object (we use ``None`` for this purpose) and checks whether the output is ``False``.
   
.. literalinclude:: ../../aeolis/tests/test_utils.py
    :language: python
    :lines: 30-45

.. note:: 
    - Classes are a way to group  and organize tests in Pytest. In the example above, the class ``TestIsIterable`` is grouping the test cases for the ``is_iterable`` function. The name of the class can be arbitrary, but it is a good practice to name it after the function that it is testing and add *Test* as prefix.
    - The ``assert`` statements are the ones that perform the actual checks. A test function can contain several ``assert`` statements. However, a best practice is to have one  or two ``assert`` statements per test function. If a test function ends up with too many ``assert`` statements, it is best to split it into two or more test functions.


Types of tests
--------------

Software tests can be of different types based on their scope of source code they test:

- **Unit tests:** test a single function or a small piece of the source code. To learn how to write unit tests for aeolis, read the section :doc:`unit testing <unit-testing>`. Unit tests for Aeolis are located in the ``aeolis/tests/`` directory.
- **Regression tests:**  a type of black-box tests where the software as a whole is tested by feeding it inputs and examining the outputs, and the internal structure of the software is rarely considered. AeoLiS currently has the following regression tests:

  - `test_netCDF_file_content.py <https://github.com/openearth/aeolis-python/blob/main/aeolis/tests/regression_tests/test_netCDF_file_content.py>`_: which tests whether a netCDF file is created upon a successful completing the modeling tasks.
  - `test_output_files_generation.py <https://github.com/openearth/aeolis-python/blob/main/aeolis/tests/regression_tests/test_output_files_generation.py>`_ which tests whether the content of the outputs (netCDF files) of running a model in the **current** version of AeoLiS are consistent with the content of outputs generated in previous versions of AeoLiS.

- **Integration tests:** test how parts of a modular software work together, or how the software behaves when it interacts with software it is expected to interact with. For example, AeoLiS has a few integrations tests that check if new versions of AeoLiS are compatible with different version of Python, starting with Python 3.9. These tests are not part of the source code, but they are executed automatically in the remote repository. To learn more about these tests, read the section on :ref:`automated testing <automated-testing>`.

.. _automated-testing:

Automated testing in AeoLiS
---------------------------

To catch any bugs introduced by new code changes, the test suite in Aeolis (unit tests + regression tests) runs automatically on every push to the main branch of the remote repository. This automated testing is enabled by `GitHub Actions <https://docs.github.com/en/actions/automating-builds-and-tests/building-and-testing-python>`_. Tests are run whenever code is pushed to the main branch (usually through a pull request), if any of the tests fail, merging changes are disabled until all tests are satisfied. The status of the test is displayed on GitHub using a green tick (pass) or a red cross (fail) next to the commit message. 

The file `.github/workflows/run-tests.yml <https://github.com/openearth/aeolis-python/blob/main/.github/workflows/run-tests.yml>`_ contains the instructions for automated testing. The tests run using the Pytest framework on a virtual machine (Ubuntu 20.04) using the GitHub-hosted runner and on Python versions 3.9-3.11. If you want to learn more about automated testing in GitHub, check `their documentation <https://docs.github.com/en/actions/automating-builds-and-tests>`_.

.. important:: 

    .. centered:: Software tests in a nutshell

    1. A software test is essentially a piece of code, such as Python modules or functions that are executed to verify whether a target software, for example, AeoLiS, is behaving as expected or not.
    #. Tests produce a pass or fail status as the output to indicate whether the target software is working as expected or not.
    #. There are different types of tests. The most common ones are unit tests, regression tests, and integration tests.
    #. In Pytest, a test script is a collection of functions that are prefixed with the word `test_`. The test functions call the functions in various modules within the aeolis source code and check whether the output is as expected. If the output is as expected, the test passes. If the output is not as expected, the test fails. This is the basic idea behind software tests. For an examplse, see the section *Example: formatting log messages* in :doc:`unit testing <unit-testing>` to learn how to write test functions. 

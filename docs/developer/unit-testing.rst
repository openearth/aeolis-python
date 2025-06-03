Unit testing
============

Unit testing is a method by which individual units of source code are tested to determine whether they behave as expected. Writing unit tests for your code is considered a best practice and is a great way to ensure that your code works as expected when you make changes to it or when other people use it. We use `Pytest <https://docs.pytest.org/en/7.1.x/index.html>`_ for writing unit test for AeoLiS. 


The first step for writing unit tests is to define test cases. A test case is a set of conditions which is used to determine whether a particular piece of functionality in your code is working as expected. For example, if you have a function which adds two numbers, you can define a test case to check whether the function returns the correct result when given two numbers, such as 5 and 3. When getting started with unit testing, *it is a good idea to write down the test cases you want to implement before writing the actual test code. This will help you to think about the functionality of your code and how it should behave in different situations.*
To define the test cases of a particular scenario in AeoLiS, use the `test scenario <https://github.com/openearth/aeolis-python/issues/new/choose>`_ template available in the GitHub issue tracker. 

The  second step is to write the actual test code. This is the code which checks whether a part of the source code is working as expected. A unit test often involves cheking that the output of a function matches one or more expected values.  If the output matches the expected result, the test case passes. If the actual result does not match the expected result, the test case fails. Below we provide a an example of how to write unit tests for AeoLiS.

Example: Formatting log messages
---------------------------------

The ``StreamFormatter`` class in the module ``model.py`` is used to customize the way log messages are formatted when printed to the console, and it is a subclass of the ``Formatter`` class from ``logging`` package. The code of the ``StreamFormatter`` class is shown below:

.. literalinclude:: ../../aeolis/model.py
   :language: python
   :pyobject: StreamFormatter
   :linenos:

The ``StreamFormatter`` has one method ``format()`` (line 4) which takes a ``record`` as input and returns a formatted string. The expected behavior of the ``format()`` method to change the format of log messages based on their level name (``INFO, ERROR, etc``). If the name *level* of the record is ``INFO``, the method shall return only the log message (lines 7-8) . In any other cases, the method shall return a single string containing the level name follwed by the log message (line 10).
To define unit test for the ``StreamFormatter`` class and it ``format()`` method, we can think of testing in terms of test scenarios and test cases, and descibe them as follows:


Test scenario
''''''''''''''''

We need to implement tests for the following scenario:

* **Description:** Test `StreamFormatter` class formats log messages based on their *level name*.
* **Code Reference:** `see repository <https://github.com/openearth/aeolis-python/blob/0596eef57f1a69946698cd0a4fdcb719411f0a1c/aeolis/model.py#L64>`_

For this scenario, we indentify the following test cases:

Test Cases
""""""""""""

1. Case: ``StreamFormatter`` is an instance of the Formatter class from the logging package.
   
   .. note:: 
      This test case is not strictly necessary, but it is a best practice to check whether the class we are testing is an instance of the class we expect it to be. This is because someone might modify the source code of this class in such a way that it is not a subclass of ``Formatter`` and the logger will not work. If that happens test must fail and we will know that we need fix the code.
   
   * **Preconditions:** None
   * **Test Steps:** 
       1. Create an instance of ``StreamFormatter`` class.
       2. Check whether the new instance is an instance of ``Formatter`` class.
   * **Test Data:** not applicable.
   * **Expected Result:** check returns ``True``.
   * **Postcondition:** None

2. Case: change the formatting of log messages other than INFO level.

   .. note:: 
      As we have described avobe, the ``format()`` method shall return the level name and the log message only when the level name is other than ``INFO``. We can check that by passing log records with different name levels to ``format()`` and compare the outputs. For example, we can pass a records with level names ``INFO`` and ``WARNING`` to ``format()`` and check if they are formatted a expected.
   
   * **Preconditions:** log records of type ``INFO`` and ``WARNING`` are created.
   * **Test Steps:** 
       1. Create an instance of ``StreamFormatter`` class.
       2. Reformat log record using the ``format()`` method.
       3. Check if log record with ``INFO`` level name is formatter as expected.
       4. Check if log record with ``WARNING`` level name is formatter as expected.
   * **Test Data:** not applicable.
   * **Expected Result:** For record with ``INFO`` level name, ``format()`` returns the log message. For recod with ``WARNING`` level name, ``format()`` return level name and log message.
   * **Postcondition:**  A ``StreamFormatter`` that will re-formatted log messeges according to level names.


Test Code
''''''''''''''''

The code for the test cases  described above can be grouped in a single test class using ``Pytest``, as follows:
 
.. literalinclude:: ../../aeolis/tests/test_model.py
   :language: python
   :pyobject: TestStreamFormatter
   :linenos:


The ``test_stream_formatter_parent()`` function test the first test case. It creates an instance of ``StreamFormatter`` and checks whether it is an instance of ``Formatter`` (lines 4-7). The ``test_stream_formatter_format()`` function test the second test case. First, it creates log records with different level names ``INFO`` and ``WARNING`` (lines 12-20). Then, it creates an instance of ``StreamFormatter`` (line 22) and re-frmats the log records using the ``format()`` method (lines 23-24). Finally, it checks whether the output of the ``format()`` method matches the expected formats (lines 25-26).

To run the tests, we can use the ``pytest`` command in the terminal as follows: ``pytest aeolis/tests/test_model.py::TestStreamFormatter``. The output should be similar to the following, if the tests pass:

.. code-block:: bash

   ==================================== test session starts ===================================
   platform linux -- Python 3.8.10, pytest-6.2.5, py-1.10.0, pluggy-1.0.0
   rootdir: /home/username/aeolis-python
   collected 2 items                                                               
   
   ... /test_model.py::TestStreamFormatter::test_stream_formatter_parent PASSED           [50%]
   ... /test_model.py::TestStreamFormatter::test_stream_formatter_info_level PASSED      [100%]
   
   ===================================== 2 passed in 0.01s ====================================


.. seealso:: 
   For learning more about testing check lesson by Code Refinery  onf `software testing <https://coderefinery.github.io/testing/>`_.
   For more information about ``Pytest`` and how to write unit tests, see the `Pytest documentation <https://docs.pytest.org/en/7.1.x/index.html>`_.

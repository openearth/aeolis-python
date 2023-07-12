Unit testing
============

Unit testing is method by which individual units of source code are tested to determine whether they are fit for use. Writing unit tests for your code is considered a best practice and is a great way to ensure that your code works as expected when you make changes to it or when other people use it. We use `Pytest <https://docs.pytest.org/en/7.1.x/index.html>`_ for writing unit test for AeoLiS. 


The first step for writing unit tests is to write test cases. A test case is a set of conditions which is used to determine whether a particular piece of functionality of your code is working as expected. For example, if you have a function which adds two numbers, you can write a test case to check whether the function returns the correct result when you pass it two numbers.

The next step is to write the actual test code. This is the code which will run the test cases and check whether the code is working as expected. The test code will run the test cases and check whether the actual result matches the expected result. If the actual result matches the expected result, the test case passes. If the actual result does not match the expected result, the test case fails.

To to define the test cases for a particular scenario, use the `test scenario <https://github.com/openearth/aeolis-python/issues/new/choose>`_ template available in the GitHub issue tracker. Below we provide a few examples.


The ``StreamFormatter`` class is used to customize the way log messages are formatted when they are printed to the console, and it is a subclass of the ``Formatter`` class in the ``logging`` package. The code of the ``StreamFormatter`` class is shown below:
=======
Example: Formatting log messages
---------------------------------


.. literalinclude:: ../../aeolis/model.py
   :language: python
   :pyobject: StreamFormatter
   :linenos:

The ``StreamFormatter`` has one method ``format()`` (line 4) which takes a ``record`` object as input and returns a formatted string. If the *level* of the record is ``INFO``, the method will return only the log message (lines 7-8) . In all other cases, the method returns the level of the messages and the message (line 10), in that order.

We can write unit test for the ``StreamFormatter`` class to ensure that the ``format()`` method behaves as expected (i.e., it works as expected). To this end, we can  describe the unit test to be implemented for this code using test scenarios and test cases. 

Test scenario
''''''''''''''''

We need to implement tests for the following scenario:


* **Description:** Test `StreamFormatter` class formats log messages based on their log level.
* **Code Reference:** `see repository <https://github.com/openearth/aeolis-python/blob/0596eef57f1a69946698cd0a4fdcb719411f0a1c/aeolis/model.py#L64>`_

For this scenario, we  indentify the following test cases:

Test Cases:
""""""""""""

1. Case: ``StreamFormatter`` is an instance of the Formatter class from the logging package.
   
   .. note:: This test case is not strictly necessary, but it is a good practice to check whether the class we are testing is an instance of the class you expect it to be. This is because some might modify the code in this class in such a way that, it is not a subclass of ``Formatter`` anymore. If this is the case, the test must fail and we will know that we need fix the code.
   

   * **Preconditions:** None
   * **Test Steps:** 
     1. Create an instance of `StreamFormatter` class.
     2. Check whether the new instance is an instance of `Formatter` class.
   * **Test Data:** not applicable.
   * **Expected Result:** check returns `True`.
   * **Postcondition:** None


2. Case: change the formatting of log messages other than INFO level.

   .. note:: 
      As we have described avobe, the ``format()`` method shall return the log level and the log message when the log level is not ``INFO``. We can test this by creating a log record with a log level other than ``INFO`` and checking whether the ``format()`` method returns a string contraining the log level and the log message.
  

   * **Preconditions:** log record of type `INFO` is created.
   * **Test Steps:** 
     1. Create an instance of `StreamFormatter` class.
     2. Create formatter log string with `format()` method
     3. Check if custom formatting is applied to log record when not `INFO` level.
   * **Test Data:** not applicable.
   * **Expected Result:** string is formatted according to the log level.
   * **Postcondition:** log message is formatted according to the log level.


Test Code
''''''''''''''''

The test code for the test cases above can be grouped in a single test class (using `Pytest`) as follows:

.. literalinclude:: ../../aeolis/tests/test_model.py
   :language: python
   :pyobject: TestStreamFormatter
   :linenos:


Explanation:
""""""""""""""

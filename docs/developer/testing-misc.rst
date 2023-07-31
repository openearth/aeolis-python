Guidelines for writing tests
----------------------------
- Start by using the 
- Write tests for the functionality that you are adding/changing
- Write tests for the functionality that should remain intact irrespective of the new code changes
- Write tests for the functionality that is affected by the changes
- Quantify the expected behavior of the function/piece of code yoou are writing the test for. For example, if you are writing a test for a function that adds two numbers, the expected behavior is that the output should be the sum of the two numbers. If you are writing a test for a function that generates a file, the expected behavior is that the file should be generated.
- Start with documenting the test design and expected behavior of the function/piece of code you are writing the test for. This will help you write the test. Use the issue tracker to document the test design and expected behavior. Use the test template to do this. Refer the issue template for the end to end test for an example.


Reviewer guidelines
-------------------

- Each pull reques to main will make the tests run automatically on the contributor's changes. As a reviewer, you should check the test results to ensure that the new changes are not breaking the existing functionality.

- Check the pass/fail status of the tests in the pull request.

- If the tests fail, the new changes are not pushed to main. The contributor should fix the issues and make a new pull request.

- If the tests fail, look into the Github actions logs to see the errors.

This page is a getting started guide on software tests for the developers of AeoLiS. It explains what software tests are, why they are important, and the different types of software tests. Read this if you are new to software tests. Skip to the next page if you are already familiar with the basic concept of software tests and want to read about the software tests developed for AeoLiS.  

Functional tests are a kind of black-box tests. They test the functionality of the model without looking at the code. The tests are based on the input and output files of the model. The input files are the same for all versions of the model. The output files are compared between the latest version of the model and the previous version of the model. If the output files across different versions of the software are the same for a given input, the test is successful. If the output files are different, the test fails.


Old script for testing introduction
-----------------------------------

Software tests are essentially a piece of code, for example, Python scripts, that are run to check whether a target software, for example AeoLiS, is working as expected or not. The tests generate a Pass or Fail status as the output upon execution to indicate this. 


- Tests are important because they help us catch bugs in the source code and ensure that the software is working as expected. 

- By running the tests as an additional set of scrips alongside the target software upon every new change to the code, we can ensure that the software is working as expected even after the changes are made. This is done using automated testing (link to section below).

- To do: Add screenshot of tests output of aeolis from the command line to show how the tests are run and what the output looks like. Also add a code snippet of the code of a test


Types of tests
--------------

- Unit tests: test a single function or a small piece of code. See the page unit testing to see the unit tests in AeoLiS.

- Integration tests: test how different pieces of code work together    

- Functional tests: Functional testing is a type of black-box testing where the software as a whole is tested by feeding it input and examining the output, and internal program structure is rarely considered. Functional testing usually describes what the system does. Functional testing does not imply that you are testing a function (method) of a module or class. Functional testing tests a slice of functionality of the whole system. 

- End-to-end tests: test the entire software


Software tests in aeolis
------------------------

aeolis currently has two types of tests:

- Unit tests: test the functions in the `aeolis` package
- Functional tests: treats the software as a black box and tests the functionality of the software as a whole by checking its output against a particular input.

Read the implementation details of the above tests in aeolis on the pages [Unit testing](unit_testing.md) and [Functional testing](functional_testing.md).
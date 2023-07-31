Reviewer guidelines
===================

It is the code maintainers' responsibility to ensure that the code is well written and well documented. This page provides guidelines for reviewing code contributions.


<Write in the same way as numpy (https://numpy.org/devdocs/dev/reviewer_guidelines.html#reviewer-guidelines)>

- Each pull reques to main will make the tests run automatically on the contributor's changes. As a reviewer, you should check the test results to ensure that the new changes are not breaking the existing functionality.

- Check if code coverage hasn't gone down. If it has, ask the contributor to add tests for the new code. 

- Check if the new code is well documented. If not, ask the contributor to add documentation. See the docstrings guide in the page quickstart.

- Check if the new code is well formatted. If not, ask the contributor to format the code. Ask the contributor to run `black` on the code.

- Add additional person as reviewer if you think they should review the code or can help review the code.
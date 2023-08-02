Automated testing in aeolis
===========================

To catch any bugs introduced by new code changes, the test suite in Aeolis (unit tests + functional tests) runs automatically on every push to the main branch. This automated testing is enabled in the aeolis repository using GitHub Actions. If the tests fail, the new changes are not pushed to main. The test outout status is displayed using a green tick (pass) or a red cross (fail) next to the commit message. 

The file .github/workflows/python-app.yml contains the instructions for the automated testing. The tests are run using the pytest framework. The tests are run on a virtual machine (Ubuntu 20.04) using the GitHub-hosted runner and on Python versions 3.8-3.11. This workflow is configured using GitHub Actions to run on every push to the main branch. 

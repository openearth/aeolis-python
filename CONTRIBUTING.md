# Contributing guidelines

Any kind of contribution to AeoLiS is welcome, from a simple comment or a question, to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). 

A contribution can be associated with the following cases:

- You have a question.
- You think you may have found a bug, including unexpected behavior.
- You want to make changes to the code base to fix a bug, make improvements, add a new functionality, or to update the documentation.

The sections below outlines the steps to make your contribution to the software for each of the aforementioned cases.

## You have a question

1. Use the search functionality [here](https://github.com/openearth/aeolis-python/issues) to see if someone already filed the same issue.
1. If your issue search did not yield any relevant results, open a new issue.
1. Apply the "Question" label. Additionally, apply other labels when relevant.

## You think you may have found a bug

1. Use the search functionality [here](https://github.com/openearth/aeolis-python/issues) to see if someone already filed the same issue.
1. If your issue search did not yield any relevant results, open a new issue and provide enough information to understand the cause and the context of the problem. Depending on the issue, you may also want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem
    - some identifying information (name and version number) for dependencies you're using
    - information about the operating system

## You want to make changes to the code base

### Announce your plan

1. (**important**) Announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue on the Github repository.
2. (**important**) Wait until a consensus is reached about your idea being a good idea.


### Set up a local development environment to work on your changes

1. If you are a part of the AeoLiS team and have write access to the aeolis-python GitHub repository, skip to the next subsection 'Develop your contribution'. All other contributors, follow the below steps.

1. Go to the [aeolis-python GitHub repository](https://github.com/openearth/aeolis-python) and click on 'Fork'. This will create a copy of the aeolis-python repository in your GitHub account. 
            
1. Clone the project to your local computer:
        
    ```bash
    git clone https://github.com/your-username/aeolis-python.git
    ```

1. Change the directory

    ```bash
    cd aeolis-python
    ```

1. Add the upstream repository

    ```bash
    git remote add upstream https://github.com/openearth/aeolis-python.git
    ```  

1. Now, git remote -v will show two remote repositories named:

        `upstream`, which refers to the aeolis-python repository <br>
        `origin`, which refers to your personal fork


### Develop your contribution

1. Create a branch off the latest commit on the AEOLIS_V2 branch to work on your feature.

    ```bash
    git checkout -b my-feature
    ```  

1. Set up a development environment on your PC by installing the package in development mode with the following command: (Consider using a virtual environment for this purpose.)

    ```bash
    pip install -e .
    ```
    
1. Set up your code editor to follow [PEP 8](https://peps.python.org/pep-0008/) (remove trailing white space, no tabs, etc.). Check code with [flake8](https://flake8.pycqa.org/en/latest/).

1. Write tests for any new lines of code you add. 

1. Include in-code documentation in form of comments and docstrings. Update the user documentation if relevant. Use the [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard) documentation style.


### Submitting your contribution


1. Push your feature branch to (your fork of) the aeolis-python GitHub repository.

1. Create a pull request, for example, following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
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

1. (**important**) Announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue on the Github repository.
1. (**important**) Wait until a consensus is reached about your idea being a good idea.
1. *Core team members with write access to the aeolis-python GitHub repository*: Clone the repository and create your own feature branch off of the latest commit on the AEOLIS_V2 branch.

    *External contributors*: Fork the repository to your own Github profile and create your own feature branch off of the latest commit on the AEOLIS_V2 branch. While working on your feature branch, make sure to stay up to date with the AEOLIS_V2 branch by pulling in changes from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. Set up a development environment on your PC by installing the package in development mode with the following command: 
    ```bash
    pip install -e .
    ```
    Consider using a virtual environment for this purpose.

1. Set up your code editor to follow [PEP 8](https://peps.python.org/pep-0008/) (remove trailing white space, no tabs, etc.). Check code with flake8.
1. Write tests for any new lines of code you add. 
1. Include in-code documentation in form of comments and docstrings. Update the user documentation if relevant. Use the [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard) documentation style.
1. Update the `CHANGELOG.md` file with the changes introduced. 
1. Push your feature branch to (your fork of) the aeolis repository on GitHub.
1. Create a pull request, for example, following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

## Acknowledgements

The contributing guidelines for AeoLiS are sourced from the [NLeSC/python-template](https://github.com/NLeSC/python-template).
Modularity
================

Our ambiton is to make the current version of AeoLiS more modular. Specially, we want the architecture of the software to consist of well defined components (modules) that integrate with each other though well defined interfaces. This will allow us to add new features more easily and to make the software more robust easier to exent.

For example, the current version of AeoLiS is organized in a single class, the `AeoLiSRunner`. This class is responsible for the execution of the model. It is also responsible for the reading of the input files, executing the solvers, and writing of the output files. This makes the code difficult to understand, extend and test, see :numref:`fig-uml-model`.

.. _fig-uml-model:

.. figure:: /images/uml-model-py.png
   :alt: AeoLiS architecture

   UML class diagram of the model.py file.


.. _GitHub repository: https://github.com/openearth/aeolis-python/issues 


Suggetions for improvement
---------------------------

We have some ideas about how to improve the modularity of the code. However, we are still working on the details of this architecture. If you have any suggestions, please let us know by opening an issue in our `GitHub repository`_ or by proposing improvements to the code via pull requests.

- **Reduce the amount to repeated code by using abstractions.** For example, some of the solvers implement similar alogorithms. We could implement an abstract class that implements the common parts of the algorithms and then extend this class to implement the specific algorithms of each solver.
- **Reduce the number of tasks the `AeoLiSRunner` and AeoLiS classes are responsible for.** This could be achieved by moving some of the code inside those clases to independent fuctions or modules. This will make the code easier to understand and test.
- **Organize code based on componets**.For example,  :numref:`fig-step-flowchart` shows a flowchart of a modeling *step* in the AeoLiSRunner. We could  use each colored-box in the flowchart to refactor the source code in separated components, each one with a well defined interface that connects to the other components. 


.. _fig-step-flowchart:

.. figure:: /images/model-step-flowchart.PNG
   :alt: AeoLiS architecture


.. seealso::
    Read about `How to write modular code` for more suggestions on how to make a software more modular.
   
   
.. _How to write modular code: https://dev.to/prxtikk/how-to-write-clean-and-modular-code-1d87


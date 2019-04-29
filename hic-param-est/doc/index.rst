Bayesian parameter estimation for relativistic heavy-ion collisions
===================================================================
.. toctree::
   :hidden:

   Home <self>

*specialized Python analysis code*

.. contents:: Outline
   :local:

Design philosophy
-----------------
As Bayesian parameter estimation becomes more common in heavy-ion physics, we need an analysis package to facilitate these projects and reduce duplication of effort.

In my experience, parameter estimation projects have many variables and "moving parts", i.e. different kinds of model calculations and experimental data, different methods for uncertainty quantification, etc.
So, I don't think it's realistic --- or necessarily even desirable --- to develop an automated analysis package.

I believe a better strategy is to create sets of modular tools for accomplishing specific projects.
The tools should be as flexible as possible without attempting to handle every conceivable corner case, and should use existing libraries, such as `scikit-learn <http://scikit-learn.org>`_, as much as possible.
Parameter estimation projects should build on previous work, adding and modifying functionality as needed, but with each project's code in its own repository.

**This means that users must have a strong understanding of Bayesian parameter estimation and must be prepared to read and modify the code.**

What is this document?
----------------------
A brief guide to https://github.com/jbernhard/hic-param-est, which is the code for my latest project applying Bayesian parameter estimation to quantify properties of the quark-gluon plasma.

I have made an effort to keep the code as generic possible so that it may be easily modified for similar projects.
Users should `fork <https://help.github.com/articles/fork-a-repo/>`_ the repository and adapt it for their needs.

Installation
------------
Python **3.5+** is required along with several common scientific Python libraries.

On most systems, the easiest way to acquire a recent version of Python and its libraries is with `Miniconda <https://conda.io/miniconda.html>`_ (be sure to use the Python 3 version).
This also has the benefit that it comes with a free version of the Intel MKL, which accelerates many linear algebra routines.
Most of the computation time in this analysis is occupied by linear algebra functions, so the MKL provides a significant speed boost.

Some Linux distributions (e.g. Arch Linux) provide the required Python packages with an optimized linear algebra library, in which case Miniconda is not necessary.

The dependencies are listed in the repository's `requirements.txt <https://github.com/jbernhard/hic-param-est/blob/master/requirements.txt>`_.
If you're using Miniconda, create and activate an environment with the dependencies::

   conda create -n hic-param-est numpy h5py scikit-learn pyyaml
   source activate hic-param-est

Install a few additional dependencies that do not have conda packages::

   pip install emcee hic

`emcee <http://dfm.io/emcee>`_ is used for MCMC sampling.
`hic <http://qcd.phy.duke.edu/hic>`_ is my library for heavy-ion collision simulation analysis, used in this project for computing flow cumulants.

Now, fork the `Github repository <https://github.com/jbernhard/hic-param-est>`_ and clone your fork.

Usage
-----
All the analysis code is in the ``src`` folder.
In the Python lexicon, each file in ``src`` is a `module <https://docs.python.org/3/tutorial/modules.html>`_ and the ``src`` folder itself is a `package <https://docs.python.org/3/tutorial/modules.html#packages>`_.

Each module is designed to be executed as a script.
Since the modules import each other, and due to the way Python intra-package references work, they must be executed using their fully-qualified module names and the ``-m`` option to the Python interpreter, e.g.::

   python -m src.expt

to run the :mod:`expt` module.
Commands should be run in the directory containing ``src``, i.e. the project root directory.

Modules
-------
.. toctree::

   modules

Stay tuned!
-----------
I plan to expand this documentation with a more complete tutorial.

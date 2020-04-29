Installation
============
TODO: Metaboverse is not yet released on conda/pypi
TODO: change package name

Conda (recommended)
-------------------

Install Miniconda, follow the steps described `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install>`_

Start the ``conda prompt``

* Windows: Open the ``Anaconda Prompt`` via the Start menu
* macOS or Linux: Open a ``Terminal``

Create a metaboverse specific ``conda`` environment.
This will install a the dependencies required to run ``metaboverse``::

    $ conda create --yes --name metaboverse -c conda-forge -c bioconda -c computational-metabolomics

.. note::

    * The installation process will take a few minutes.
    * Feel free to use a different name for the Conda environment

    You can use the following command to remove a conda environment::

        $ conda env remove -y --name metaboverse

    This is only required if something has gone wrong in the previous step.

Activate the ``metaboverse`` environment::

    $ conda activate metaboverse

To test your ``metaboverse`` installation, in your Conda Prompt, run the command::

    $ metaboverse --help

or::

    $ python
    import metaboverse

Close and deactivate the ``metaboverse`` environment when you’re done::

    $ conda deactivate


PyPi
----

Install the current release of ``metaboverse`` with ``pip``::

    $ pip install metaboverse

.. note::

    * The installation process will take a few minutes.

To upgrade to a newer release use the ``--upgrade`` flag::

    $ pip install --upgrade metaboverse

If you do not have permission to install software systemwide, you can
install into your user directory using the ``--user`` flag::

    $ pip install --user metaboverse

Alternatively, you can manually download ``metaboverse`` from
`GitHub <https://github.com/computational-metabolomics/metaboverse/releases>`_  or
`PyPI <https://pypi.python.org/pypi/metaboverse>`_.
To install one of these versions, unpack it and run the following from the
top-level source directory using the Terminal::

    $ pip install .

Testing
-------
*Metaboverse* uses the Python ``unittest`` testing package.
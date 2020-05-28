Installation
============
TODO: MetaboBlend is not yet released on conda/pypi
TODO: change package name

Conda (recommended)
-------------------

Install Miniconda, follow the steps described `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install>`_

Start the ``conda prompt``

* Windows: Open the ``Anaconda Prompt`` via the Start menu
* macOS or Linux: Open a ``Terminal``

Create a metaboblend specific ``conda`` environment.
This will install a the dependencies required to run ``metaboblend``::

    $ conda create --yes --name metaboblend -c conda-forge -c bioconda -c computational-metabolomics

.. note::

    * The installation process will take a few minutes.
    * Feel free to use a different name for the Conda environment

    You can use the following command to remove a conda environment::

        $ conda env remove -y --name metaboblend

    This is only required if something has gone wrong in the previous step.

Activate the ``metaboblend`` environment::

    $ conda activate metaboblend

To test your ``metaboblend`` installation, in your Conda Prompt, run the command::

    $ metaboblend --help

or::

    $ python
    import metaboblend

Close and deactivate the ``metaboblend`` environment when you’re done::

    $ conda deactivate


PyPi
----

Install the current release of ``metaboblend`` with ``pip``::

    $ pip install metaboblend

.. note::

    * The installation process will take a few minutes.

To upgrade to a newer release use the ``--upgrade`` flag::

    $ pip install --upgrade metaboblend

If you do not have permission to install software systemwide, you can
install into your user directory using the ``--user`` flag::

    $ pip install --user metaboblend

Alternatively, you can manually download ``metaboblend`` from
`GitHub <https://github.com/computational-metabolomics/metaboblend/releases>`_  or
`PyPI <https://pypi.python.org/pypi/metaboblend>`_.
To install one of these versions, unpack it and run the following from the
top-level source directory using the Terminal::

    $ pip install .

Testing
-------
*MetaboBlend* uses the Python ``unittest`` testing package.
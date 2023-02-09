Installing pyroe
================

The ``pyroe`` package can be accessed from its `github repository <https://github.com/COMBINE-lab/pyroe>`_, installed via `pip <https://pip.pypa.io/en/stable/>`_. To install the ``pyroe`` package via ``pip`` use the command:

.. code:: bash

  pip install pyroe


To make use of the ``load_fry`` function (which, itself, installs `scanpy <https://scanpy.readthedocs.io/en/stable/>`_):

.. code:: bash

  pip install pyroe


Alternatively, ``pyroe`` can be installed via ``bioconda``, which will automatically install the variant of the package including ``load_fry``, and will
also install ``bedtools`` to enable faster construction of the *splici* reference.  This installation can be performed with the following shell command:

.. code:: bash

  conda install pyroe -c bioconda


with the appropriate bioconda channel in the conda channel list.

Deploying ``lcdb-wf`` and staying up-to-date
============================================
The repository comes with lots of infrastructure for testing that is not
necessarily needed in practice when using lcdb-wf for a project. To get
a simplified version:

.. code-block:: bash

   python deploy.py --flavor rnaseq project-dir


The script will use ``rsync`` to copy over files to `project-dir`, excluding
various test files and excluding any files that may have been created in the
process of testing. For "flavor", choose ``chipseq``, ``rnaseq``,
``colocalization``, or ``full`` to get everything. 

This script also writes a file in the destination called
``.lcdb-wf-deployment.json`` which stores details about what commit was used to
deploy and the timestamp. This can come in handy later when comparing
a deployed directory with the main repository to decide whether to update.

Updating
--------
The most straightforward approach to updating is to use a diff tool like `meld
<http://meldmerge.org>`_ to visually compare differences between a deployed
project and a freshly-cloned version of lcdb-wf:

.. code-block:: bash

   git clone https://github.com/lcdb/lcdb-wf.git comparison-directory
   meld project-dir comparison-directory

This way you can pick and choose which updates are relevant without having to
resort to arcane git commands and difficult merges.

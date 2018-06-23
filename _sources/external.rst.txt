.. _external:

"External" workflow
-------------------
Often we want to compare new data with existing published data. We have found
that in practice, having a separate workflow to handle downloading and
reformatting and various conversion tasks helps with organization.

The test workflow is a working example that:

- downloads ChIP-seq data from modENCODE in an older fly genome organism (dm3)
- downloads the chainfile for liftover
- fixes the formatting of the downloaded files so they can be lifted over
- lifts over the files to the newer dm6 assembly.

The file is intended to be heavily edited for the particular experiment; it is
here mostly as a placeholder and to be used as a template for integrative
downstream work.  It can then be incorporated into the ``figures`` workflow (see
:ref:`figures`) to integrate the analysis with other output.

.. image:: external.png

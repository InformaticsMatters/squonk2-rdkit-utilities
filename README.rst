Informatics Matters RDKit Utilities
====================================

.. image:: https://badge.fury.io/py/im-rdkit-utilities.svg
   :target: https://badge.fury.io/py/im-rdkit-utilities
   :alt: PyPI package (latest)

.. image:: https://github.com/InformaticsMatters/squonk2-rdkit-utilities/actions/workflows/build.yaml/badge.svg
   :target: https://github.com/InformaticsMatters/squonk2-rdkit-utilities/actions/workflows/build.yaml
   :alt: Build

.. image:: https://github.com/InformaticsMatters/squonk2-rdkit-utilities/actions/workflows/publish.yaml/badge.svg
   :target: https://github.com/InformaticsMatters/squonk2-rdkit-utilities/actions/workflows/publish.yaml
   :alt: Publish

A Python package of RDKit-specific helpers shared by **Squonk2 Data Manager
Jobs**: molecule readers/writers over SDF and delimited-SMILES text formats,
fragment selection, and a handful of small molecule-inspection helpers.

This consolidates the ``rdkit_utils.py`` module that had been copy-pasted,
and had begun to diverge, across several Job repositories
(``squonk2-desc-rdkit``, ``squonk2-desc-mordred``, ``squonk2-jaqpot`` and
``virtual-screening``) into a single, tested, authoritative source.

Installation (Python)
======================

The package is published on `PyPI`_ and can be installed from there::

    pip install im-rdkit-utilities

Once installed, import it as ``rdkit_utils`` (matching the module name it
replaces, so existing call sites need no changes beyond the import path)::

    >>> import rdkit_utils
    >>> reader = rdkit_utils.create_reader('molecules.smi', delimiter='\t')

Public surface
==============

- ``create_reader()`` / ``create_writer()`` — construct a reader/writer for
  a .sdf, .sdf.gz or delimited-SMILES file, based on its extension.
- ``SdfReader`` / ``SdfWriter`` / ``SmilesReader`` / ``SmilesWriter`` — the
  underlying reader/writer implementations.
- ``generate_headers()`` — build output headers for a tab/comma separated
  file, given the ID column configuration.
- ``fragment()`` — pick the largest fragment of a (typically salted)
  molecule, by heavy-atom count or molecular weight.
- ``fragmentAndFingerprint()`` — fragment a stream of molecules and
  fingerprint the result.
- ``get_num_chiral_centers()`` / ``get_num_sp3_centres()`` — small
  molecule-inspection helpers.
- ``check_molecules_are_3d()`` — check whether the molecules in a .sdf file
  have 3D conformers.
- ``rdk_read_single_mol()`` / ``rdk_read_mols()`` / ``rdk_read_molecule_files()``
  / ``rdk_merge_mols()`` / ``rdk_mol_supplier()`` / ``sdf_record_gen()`` —
  molecule-file reading helpers.
- ``updateChargeFlagInAtomBlock()`` — adds the legacy charge-flag encoding
  to a full CTAB molblock (counts line included), for tools such as rDock
  that only understand the old syntax. This is **not** the same as
  ``dm_job_utilities.utils.update_charge_flag_in_atom_block()``, which
  operates on just the atom-block portion of a molblock (one line further
  in) — the two take different input shapes and are not interchangeable.
  This RDKit-oriented variant lives here because it's used directly
  alongside the readers/writers above; the job-utilities function remains
  the natural home for pure string manipulation with no RDKit dependency.

Command-line helpers
====================

- ``add_common_molecule_io_args(parser)`` — adds the "Input/output options"
  argument group shared by the molecule processing Jobs (``-i/--infile``,
  ``-o/--outfile``, ``-d/--delimiter``, ``--id-column``, ``--mol-column``,
  ``--read-header``, ``--write-header``, ``--read-records``,
  ``-k/--omit-fields``, and optionally ``--y-column``), and returns the group
  so more options can be added to it. Keyword arguments ``output_default``,
  ``output_required`` and ``include_y_column`` cover the variations between
  Jobs.
- ``str_or_int()`` — an argparse ``type`` for a column specifier given either
  as a zero-based index or as a field name.

The group's namespace feeds ``create_reader()`` directly::

    >>> import argparse, rdkit_utils
    >>> from dm_job_utilities.utils import read_delimiter
    >>> parser = argparse.ArgumentParser()
    >>> _ = rdkit_utils.add_common_molecule_io_args(parser)
    >>> args = parser.parse_args(['-i', 'molecules.smi', '-d', 'tab'])
    >>> reader = rdkit_utils.create_reader(
    ...     args.input,
    ...     delimiter=read_delimiter(args.delimiter),
    ...     read_header=args.read_header,
    ...     id_column=args.id_column,
    ...     mol_column=args.mol_column,
    ...     read_records=args.read_records,
    ... )

``--infile`` and ``--outfile`` are the canonical long spellings. The parsed
values are always available as ``args.input`` and ``args.output``.

``--input`` and ``--output`` are **deprecated aliases**, retained only so that
adopting the helper did not break existing Job manifests. They are due to be
removed — see `issue #5
<https://github.com/InformaticsMatters/squonk2-rdkit-utilities/issues/5>`_. Job
Definitions should invoke ``-i``/``-o``: they are canonical, unaffected by that
removal, and work against container images built before the Job adopted the
helper.

Defaults follow the ``create_reader()`` / ``create_writer()`` signatures rather
than the values the Job scripts hand-type today — in particular ``--mol-column``
defaults to ``None`` so that ``SmilesReader`` can infer it, where most Jobs
currently pass an explicit ``0``.

Progress and cost reporting (``--interval``, ``ProgressReporter``) is
deliberately **not** here — that is a Data Manager logging concern and lives in
`im-data-manager-job-utilities`_.

.. _PyPI: https://pypi.org/project/im-rdkit-utilities
.. _im-data-manager-job-utilities: https://pypi.org/project/im-data-manager-job-utilities

Get in touch
============

- Report bugs, suggest features or view the source code `on GitHub`_.

.. _on GitHub: https://github.com/informaticsmatters/squonk2-rdkit-utilities

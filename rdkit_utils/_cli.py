"""Command-line helpers for Squonk2 Data Manager Jobs that read and write
molecule files.

The "Input/output options" argparse group defined here is currently re-typed,
with small divergences, across the Job scripts that use this package's readers
and writers. Its options correspond directly to the ``create_reader()`` and
``create_writer()`` parameters in ``_core``, which is why the builder lives
beside them rather than in ``dm_job_utilities``.

Progress and cost reporting (``--interval``, ``ProgressReporter``) is
deliberately *not* here. That is a Data Manager logging concern with nothing
RDKit about it, and lives in ``dm_job_utilities``.
"""

import argparse
from typing import Optional, Union


def str_or_int(value: str) -> Union[str, int]:
    """An argparse ``type`` for a column specifier that may be given either as
    a zero-based integer index or as a field name.

    Which form is meaningful depends on the file being read: ``SmilesReader``
    indexes columns by position, whereas ``SdfReader`` matches an SDF property
    by name. Returns an ``int`` when the value looks like one, otherwise the
    original string.
    """
    try:
        return int(value)
    except ValueError:
        return value


def add_common_molecule_io_args(
    parser: argparse.ArgumentParser,
    *,
    output_default: Optional[str] = None,
    output_required: bool = False,
    include_y_column: bool = False,
) -> argparse._ArgumentGroup:  # pylint: disable=protected-access
    """Add the "Input/output options" argument group shared by the molecule
    processing Jobs, and return the group so a caller can add more options to
    it.

    Adds ``-i/--infile``, ``-o/--outfile``, ``-d/--delimiter``,
    ``--id-column``, ``--mol-column``, ``--read-header``, ``--write-header``,
    ``--read-records`` and ``-k/--omit-fields``; ``--y-column`` is added when
    *include_y_column* is set.

    ``--infile`` and ``--outfile`` are the canonical spellings. ``--input`` and
    ``--output`` are accepted as aliases so that adopting this helper does not
    break existing Job manifests, and the parsed values remain available as
    ``args.input`` and ``args.output``.

    The defaults deliberately match the ``create_reader()`` / ``create_writer()``
    signatures rather than the values the Job scripts hand-type today:

    - ``--mol-column`` defaults to ``None``, not ``0``. ``SmilesReader`` infers
      the molecule column when it is not given (``1`` if the ID column is
      ``0``, otherwise ``0``); passing an explicit ``0`` defeats that.
    - ``--delimiter`` defaults to ``None``, which ``SmilesReader`` treats as
      whitespace. Pass the value through
      ``dm_job_utilities.utils.read_delimiter()`` first to resolve the symbolic
      names ``tab``, ``space``, ``comma`` and ``pipe``.
    - ``--read-records`` defaults to ``100``, matching ``create_reader()``.

    :param output_default: default for ``--outfile`` when the Job writes to a
        conventional file name unless told otherwise.
    :param output_required: require ``--outfile``. Off by default; several Jobs
        need it on.
    :param include_y_column: add ``--y-column``, for the Jobs that model
        against a response variable.
    """
    group = parser.add_argument_group("Input/output options")
    group.add_argument(
        "-i",
        "--infile",
        "--input",
        dest="input",
        required=True,
        help="Input file (.smi, .txt, .sdf or .sdf.gz)",
    )
    group.add_argument(
        "-o",
        "--outfile",
        "--output",
        dest="output",
        default=output_default,
        required=output_required,
        help="Output file (.smi, .txt or .sdf)",
    )
    # to pass tab as the delimiter specify it as $'\t' or use one of the
    # symbolic names 'comma', 'tab', 'space' or 'pipe'
    group.add_argument("-d", "--delimiter", help="Delimiter when using SMILES")
    group.add_argument(
        "--id-column",
        type=str_or_int,
        help="Column for the molecule ID"
        " (zero-based index for .smi, field name for .sdf)",
    )
    group.add_argument(
        "--mol-column",
        type=str_or_int,
        help="Column for the molecule when using delimited text formats"
        " (zero-based index)",
    )
    if include_y_column:
        group.add_argument(
            "--y-column",
            type=str_or_int,
            help="Column for the Y variable"
            " (zero-based index for .smi, field name for .sdf)",
        )
    group.add_argument(
        "--read-header",
        action="store_true",
        help="Read a header line with the field names when reading .smi or .txt",
    )
    group.add_argument(
        "--write-header",
        action="store_true",
        help="Write a header line when writing .smi or .txt",
    )
    group.add_argument(
        "--read-records",
        type=int,
        default=100,
        help="Read this many records to determine the fields that are present",
    )
    group.add_argument(
        "-k",
        "--omit-fields",
        action="store_true",
        help="Don't include fields from the input in the output",
    )
    return group

"""RDKit utilities for Squonk2 Data Manager Jobs.

Everything in ``_core`` is re-exported here so callers can keep using the
flat ``import rdkit_utils`` / ``rdkit_utils.create_reader(...)`` style that
predates this package.
"""

from rdkit_utils._cli import add_common_molecule_io_args, str_or_int
from rdkit_utils._core import (
    ID_COL_NAME,
    SMILES_COL_NAME,
    SdfReader,
    SdfWriter,
    SmilesReader,
    SmilesWriter,
    check_molecules_are_3d,
    create_reader,
    create_writer,
    fragment,
    fragmentAndFingerprint,
    generate_headers,
    get_num_chiral_centers,
    get_num_sp3_centres,
    rdk_merge_mols,
    rdk_mol_supplier,
    rdk_read_molecule_files,
    rdk_read_mols,
    rdk_read_single_mol,
    sdf_record_gen,
    updateChargeFlagInAtomBlock,
)

__all__ = [
    "ID_COL_NAME",
    "SMILES_COL_NAME",
    "SdfReader",
    "SdfWriter",
    "SmilesReader",
    "SmilesWriter",
    "add_common_molecule_io_args",
    "check_molecules_are_3d",
    "create_reader",
    "create_writer",
    "fragment",
    "fragmentAndFingerprint",
    "generate_headers",
    "get_num_chiral_centers",
    "get_num_sp3_centres",
    "rdk_merge_mols",
    "rdk_mol_supplier",
    "rdk_read_molecule_files",
    "rdk_read_mols",
    "rdk_read_single_mol",
    "sdf_record_gen",
    "str_or_int",
    "updateChargeFlagInAtomBlock",
]

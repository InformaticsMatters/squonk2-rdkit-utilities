"""RDKit-specific utilities for Squonk2 Data Manager Jobs.

Molecule readers and writers over SDF and delimited-SMILES text formats,
fragment selection, and a handful of small molecule-inspection helpers.

Consolidated from the four near-identical ``rdkit_utils.py`` copies that had
accumulated across ``squonk2-desc-rdkit``, ``squonk2-desc-mordred``,
``squonk2-jaqpot`` and ``virtual-screening``, using the desc-rdkit/
desc-mordred pair (identical, and the most complete of the four) as the
canonical base.
"""

import csv
import gzip
from typing import Any, List, Optional

from dm_job_utilities.dm_log import DmLog
from dm_job_utilities.utils import log
from rdkit import Chem
from rdkit.Chem import Descriptors

ID_COL_NAME = "ID"
SMILES_COL_NAME = "SMILES"


class SdfWriter:
    """Writes molecules, with calculated properties, to a .sdf/.sdf.gz file."""

    def __init__(self, outfile: str, prop_names: List[str]):
        self.gzip: Optional[Any] = None
        if outfile.endswith(".gz"):
            self.gzip = gzip.open(outfile, "wt")
            self.writer = Chem.SDWriter(self.gzip)
        else:
            self.writer = Chem.SDWriter(outfile)
        self.prop_names = prop_names

    def write(
        self, smi, mol, mol_id, existing_props, new_props, smiles_prop_name=None
    ):
        """Write a molecule and its calculated properties to the SDF."""
        del existing_props
        if not mol:
            mol = Chem.MolFromSmiles(smi)
        if mol_id is not None:
            mol.SetProp("_Name", mol_id)
        if smiles_prop_name is not None:
            mol.SetProp(smiles_prop_name, smi)

        for index, prop_name in enumerate(self.prop_names):
            value = new_props[index]
            if prop_name is not None:
                mol.SetProp(prop_name, str(value))

        self.writer.write(mol)

    def write_header(self, values):
        """SDF has no header line to write; this is a documented no-op."""
        del values
        log("INFO: asked to write header for an SDF. No action will be taken.")

    def close(self):
        """Close the underlying SDF (and, if used, gzip) writer."""
        self.writer.close()
        if self.gzip:
            self.gzip.close()


class SmilesWriter:
    """Writes molecules to a delimited SMILES (.smi/.txt) file."""

    def __init__(self, outfile, sep, extra_field_names, id_column=None, mol_column=0):
        # The handle is held open for the writer's lifetime (write() is
        # called once per molecule); 'with' doesn't apply here.
        self.writer = open(  # pylint: disable=consider-using-with
            outfile, "w", encoding="utf-8"
        )
        self.sep = " " if sep is None else sep
        self.extra_field_names = extra_field_names

        self.id_column = None if id_column is None else int(id_column)
        # Mirror SmilesReader's inference, so a writer built from the same
        # arguments as its reader lays the columns out the same way. write()
        # orders the ID and SMILES columns by comparing the two, so this must
        # not be left as None.
        if mol_column is None:
            self.mol_column = 1 if self.id_column == 0 else 0
        else:
            self.mol_column = int(mol_column)

    def write_header(self, values):
        """Write a header line built from the given column names."""
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def write(
        self, smi, mol, mol_id, existing_props, new_props, smiles_prop_name=None
    ):
        """Write a molecule (as SMILES), its ID and properties as a row."""
        del mol, smiles_prop_name
        if isinstance(mol_id, str) and self.sep in mol_id:
            mol_id = '"' + mol_id + '"'

        if self.id_column is not None:
            if self.id_column < self.mol_column:
                values = [mol_id, smi]
            else:
                values = [smi, mol_id]
        else:
            values = [smi]

        for prop in existing_props:
            if prop is not None:
                if isinstance(prop, str) and self.sep in prop:
                    values.append('"' + prop + '"')
                else:
                    values.append(prop)
            else:
                values.append("")

        for prop in new_props:
            values.append("" if prop is None else str(prop))
        line = self.sep.join(values)
        self.writer.write(line + "\n")

    def close(self):
        """Close the underlying writer."""
        self.writer.close()


class SdfReader:
    """Reads molecules from a .sdf/.sdf.gz file."""

    def __init__(self, filename, id_col, recs_to_read):

        self.field_names: List[str] = []
        if id_col == "_Name":
            self.field_names.append(ID_COL_NAME)
        # read a number of records to determine the field names
        if recs_to_read:
            self._collect_field_names(filename, recs_to_read)

        # now create the real reader
        self.reader = self.create_reader(filename)
        self.id_col = id_col

    def _collect_field_names(self, filename, recs_to_read):
        """Scan up to recs_to_read records to discover the SDF's property
        names, so they can be used as output field names later.
        """
        reader = self.create_reader(filename)
        for _ in range(recs_to_read):
            if reader.atEnd():
                break
            try:
                mol = next(reader)
            except StopIteration:
                break
            if not mol:
                continue
            for name in mol.GetPropNames():
                if name not in self.field_names:
                    self.field_names.append(name)

    def get_mol_field_name(self):
        """SDF readers have no distinct 'molecule' column name."""
        return None

    @staticmethod
    def create_reader(filename):
        """Create a ForwardSDMolSupplier for a .sdf or .sdf.gz file."""
        if filename.endswith(".gz"):
            return Chem.ForwardSDMolSupplier(gzip.open(filename))
        return Chem.ForwardSDMolSupplier(filename)

    def read(self):
        """Read the next molecule, returning (mol, smiles, mol_id, props),
        or None once the file is exhausted.
        """
        try:
            props = []
            mol = next(self.reader)
            if mol is None:
                return None, None, None, None
            smi = Chem.MolToSmiles(mol)
            if self.id_col:
                mol_id = mol.GetProp(self.id_col) if mol.HasProp(self.id_col) else None
            else:
                mol_id = None
            if self.id_col == "_Name":
                props.append(mol_id)

            for name in self.field_names:
                props.append(mol.GetProp(name) if mol.HasProp(name) else None)
            return mol, smi, mol_id, props

        except StopIteration:
            return None

    def get_extra_field_names(self):
        """Return the SDF property names discovered when the reader was
        created (or supplied by the caller via the id_col field).
        """
        return self.field_names

    def close(self):
        """No-op: ForwardSDMolSupplier has no explicit close."""


class SmilesReader:
    """Reads molecules from a delimited SMILES (.smi/.txt) file."""

    def __init__(
        self, filename, read_header, delimiter, id_column, mol_column, recs_to_read
    ):
        self.delimiter = delimiter if delimiter is not None else " "
        if mol_column is None:
            self.mol_column = 1 if id_column == 0 else 0
        else:
            self.mol_column = mol_column

        self.id_column = None if id_column is None else int(id_column)

        tmp_reader, file_reader = self.create_readers(filename)
        # read header line
        if read_header:
            tokens = next(tmp_reader)
            self.field_names = [token.strip() for token in tokens]
        else:
            self.field_names = [None] * max(
                1,
                self.mol_column + 1,
                0 if self.id_column is None else self.id_column + 1,
            )
            self.field_names[self.mol_column] = SMILES_COL_NAME
            if self.id_column is not None:
                self.field_names[self.id_column] = ID_COL_NAME

        max_num_tokens = 0
        for _ in range(0, recs_to_read):
            line = next(tmp_reader, None)
            if not line:
                break
            max_num_tokens = max(max_num_tokens, len(line))

        if max_num_tokens > len(self.field_names):
            for index in range(len(self.field_names), max_num_tokens):
                self.field_names.append("field" + str(index + 1))

        file_reader.close()

        # now create the real reader and discard the header
        self.reader, self.file = self.create_readers(filename)
        if read_header:
            next(self.reader)

    def create_readers(self, filename):
        """Open filename (transparently gunzipping) and wrap it in a csv
        reader, returning (csv_reader, file_handle).
        """
        # The handle is held open for the reader's lifetime (read() is
        # called once per record); 'with' doesn't apply here.
        # pylint: disable=consider-using-with
        if filename.endswith(".gz"):
            handle = gzip.open(filename, "rt", encoding="utf-8")
        else:
            handle = open(filename, "rt", encoding="utf-8")
        return csv.reader(handle, delimiter=self.delimiter), handle

    def get_mol_field_name(self):
        """Return the field name of the molecule (SMILES) column."""
        return self.field_names[self.mol_column] if self.field_names else None

    def read(self):
        """Read the next record, returning (mol, smiles, mol_id, props),
        or None once the file is exhausted.
        """
        tokens = next(self.reader, None)
        if not tokens:
            return None
        smi = tokens[self.mol_column]
        mol_id = tokens[self.id_column] if self.id_column is not None else None

        mol = Chem.MolFromSmiles(smi)
        props = []

        for index, token in enumerate(tokens):
            token = token.strip()
            if index not in (self.mol_column, self.id_column):
                props.append(token)
                if mol:
                    if self.field_names and len(self.field_names) > index:
                        mol.SetProp(self.field_names[index], token)
                    else:
                        mol.SetProp("field" + str(index), token)

        return mol, smi, mol_id, props

    def get_extra_field_names(self):
        """Return the input's field names, excluding the molecule and ID
        columns.
        """
        if not self.field_names:
            return []
        return [
            name
            for index, name in enumerate(self.field_names)
            if index not in (0, self.id_column)
        ]

    def close(self):
        """Close the underlying file handle."""
        self.file.close()


def generate_headers(
    id_col_type,
    id_col_value,
    mol_field_name: Optional[str],
    field_names: List[str],
    calc_prop_names: List[str],
    omit_fields: bool,
):
    """Generate the headers for when writing a tab or comma separated file.

    :param id_col_type: The type of ID column that was specified.
        -1 String for the SDF field name
        +1 int for column index for TAB
        0 for None (no ID col)
    :param id_col_value: The value specified for the ID column
        -1: the string specified
        +1: the value specified as an int
        0: None
    :param mol_field_name: The name of the mol column, or None for SDF
    :param field_names: The names of the fields in the input to be reused
    :param calc_prop_names: The names of the new fields to be added
    :param omit_fields: Do not add the input fields (except for mol and ID)
    :return: List of the header names
    """
    headers = []
    mol_header = mol_field_name if mol_field_name else SMILES_COL_NAME

    if id_col_type == 0:  # id_col was None
        if omit_fields or mol_header not in field_names:
            headers.append(mol_header)
    elif id_col_type == 1:  # was an int so CSV
        if omit_fields or mol_header not in field_names:
            headers.append(mol_header)
        if omit_fields:
            headers.append(field_names[id_col_value])
    else:  # id_col was string, so SDF field name or _Name
        headers.append(mol_header)
        id_name = ID_COL_NAME if id_col_value == "_Name" else id_col_value
        if id_name not in field_names or omit_fields:
            headers.append(id_name)

    if not omit_fields:
        headers.extend(field_names)

    headers.extend(calc_prop_names)
    return headers


def create_reader(
    filename,
    filetype=None,
    id_column=None,
    mol_column=None,
    read_records=100,
    read_header=False,
    delimiter="\t",
):
    """Create a reader for the given file, choosing SdfReader or SmilesReader
    based on filetype (or the filename's extension when filetype is None).
    """
    if filetype is None:
        if filename.endswith((".sdf", ".sdf.gz", ".sd", ".sd.gz")):
            filetype = "sdf"
        else:
            filetype = "smi"

    if filetype == "sdf":
        return SdfReader(filename, id_column, read_records)
    if filetype == "smi":
        return SmilesReader(
            filename, read_header, delimiter, id_column, mol_column, read_records
        )
    raise ValueError("Unexpected file type", filetype)


def create_writer(
    outfile,
    delimiter="\t",
    extra_field_names=None,
    calc_prop_names=None,
    id_column=None,
    mol_column=0,
):
    """Create a writer for the given file, choosing SdfWriter or SmilesWriter
    based on the outfile's extension.
    """
    if not extra_field_names:
        extra_field_names = []
    if not calc_prop_names:
        calc_prop_names = []
    if outfile.endswith(".sdf") or outfile.endswith("sd"):
        return SdfWriter(outfile, calc_prop_names)
    return SmilesWriter(
        outfile, delimiter, extra_field_names, id_column=id_column, mol_column=mol_column
    )


def updateChargeFlagInAtomBlock(mb):  # pylint: disable=invalid-name
    """Add data for the charges to a full CTAB molblock (counts line included).

    This data is deprecated and should be specified using "M  CHG" lines,
    but some old software such as rDock only handles the old syntax. RDKit
    only supports the new syntax so this function adds the old syntax back
    in alongside it. Based on work by Jose Manuel Gally, see
    https://sourceforge.net/p/rdkit/mailman/message/36425493/

    Note this operates on a *molblock* (starting with the title/counts
    lines), whereas ``dm_job_utilities.utils.update_charge_flag_in_atom_block``
    operates on just the *atom block* portion (one line further in) — the
    two are not interchangeable. This variant is kept here, alongside the
    RDKit-oriented readers/writers that produce full molblocks; the
    job-utilities function remains the natural home for pure atom-block
    string manipulation with no RDKit dependency.
    """
    formatter = "{:>10s}" * 3 + " {:>2}" + "{:>3s}" * 12
    chgs = []  # list of (atom index, charge) tuples
    lines = mb.split("\n")
    ctab = lines[3]
    atom_count = int(ctab.split()[0])
    # parse the molblock line by line
    for line in lines:
        # look for the M CHG property
        if line[0:6] == "M  CHG":
            # "M  CHG X" is not needed for parsing; the values we want follow
            records = line.split()[3:]
            for index in range(0, len(records), 2):
                idx = records[index]
                chg = records[index + 1]
                chgs.append((int(idx), int(chg)))
            break

    # sort by idx so the molblock only needs parsing once more
    chgs = sorted(chgs, key=lambda x: x[0])

    # now attribute each charge to its atom line
    for chg in chgs:
        index = 4
        while index < 4 + atom_count:
            if index - 3 == chg[0]:
                fields = lines[index].split()
                x = fields[0]  # pylint: disable=invalid-name
                y = fields[1]  # pylint: disable=invalid-name
                z = fields[2]  # pylint: disable=invalid-name
                symb = fields[3]
                mass_diff = fields[4]
                charge = fields[5]
                sp = fields[6]  # pylint: disable=invalid-name
                hc = fields[7]  # pylint: disable=invalid-name
                scb = fields[8]
                v = fields[9]  # pylint: disable=invalid-name
                hd = fields[10]  # pylint: disable=invalid-name
                nu1 = fields[11]
                nu2 = fields[12]
                aamn = fields[13]
                irf = fields[14]
                ecf = fields[15]
                if chg[1] == -1:
                    charge = "5"
                elif chg[1] == -2:
                    charge = "6"
                elif chg[1] == -3:
                    charge = "7"
                elif chg[1] == 1:
                    charge = "3"
                elif chg[1] == 2:
                    charge = "2"
                elif chg[1] == 3:
                    charge = "1"
                else:
                    log(
                        "ERROR! "
                        + str(lines[0])
                        + "unknown charge flag: "
                        + str(chg[1])
                    )
                    charge = "0"
                lines[index] = formatter.format(
                    x, y, z, symb, mass_diff, charge, sp, hc, scb, v, hd,
                    nu1, nu2, aamn, irf, ecf,
                )
            index += 1
    return "\n".join(lines)


def get_num_chiral_centers(mol):
    """Return (number of chiral centres, number of those left undefined)."""
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = sum(1 for center in chiral_centers if center[1] == "?")
    return len(chiral_centers), undef_cc


def get_num_sp3_centres(mol):
    """Return the number of sp3-hybridised atoms in the molecule."""
    return sum(
        atom.GetHybridization() == Chem.HybridizationType.SP3 for atom in mol.GetAtoms()
    )


def sdf_record_gen(hnd):
    """A generator for text records from a SD file."""
    mol_text_tmp = ""
    while True:
        line = hnd.readline()
        if not line:
            return
        mol_text_tmp += line.decode("utf-8")
        if line.decode("utf-8").startswith("$$$$"):
            yield mol_text_tmp
            mol_text_tmp = ""


def rdk_read_single_mol(filename):
    """Read a single molecule. Can be a .mol, .sdf or .sdf.gz file.
    If a .sdf, the first molecule is returned.
    """
    if filename.endswith(".mol"):
        return Chem.MolFromMolFile(filename)
    if filename.endswith(".sdf"):
        return next(Chem.SDMolSupplier(filename))
    if filename.endswith(".sdf.gz"):
        with gzip.open(filename, "rb") as gz:
            return next(Chem.ForwardSDMolSupplier(gz))
    raise ValueError(
        "Unsupported file type. Must be .mol .sdf or .sdf.gz. Found " + filename
    )


def rdk_read_mols(filename):
    """Read molecules from a single file (.mol, .sdf or .sdf.gz)."""
    if filename.endswith(".mol"):
        return [Chem.MolFromMolFile(filename)]
    return list(rdk_mol_supplier(filename))


def rdk_read_molecule_files(inputs):
    """Read input molecules.

    A list of inputs is specified, each one containing either a single
    filename or a comma separated list of filenames (no spaces). The files
    can either be .mol files with a single molecule or .sdf files with
    multiple molecules. The molecules in all the files specified are read
    and returned as a list of molecules.

    Examples::

        ['mols.sdf']                          - a single SDF, one or more mols
        ['mol1.mol', 'mol2.mol', 'mol3.mol']  - 3 molfiles, separate elements
        ['mol1.mol,mol2.mol,mol3.mol']        - 3 molfiles, one CSV element
        ['mol1.mol,mol2.mol', 'mols.sdf']     - 3 molfiles plus 1 SDF

    :param inputs: Input filenames
    :return: List of molecules that have been read
    """
    mols = []
    for input_names in inputs:
        for token in input_names.split(","):
            if token.endswith(".mol"):
                mol = Chem.MolFromMolFile(token)
                if mol:
                    mols.append(mol)
                else:
                    DmLog.emit_event("WARNING: could not process", token)
            else:
                for index, mol in enumerate(rdk_read_mols(token)):
                    if mol:
                        mols.append(mol)
                    else:
                        DmLog.emit_event(
                            "WARNING: could not process molecule", index, "from", token
                        )
    return mols


def rdk_merge_mols(inputs):
    """Merge multiple molecules into a single molecule.

    :param inputs: Input filenames, in the format accepted by
        ``rdk_read_molecule_files``.
    :return: Tuple of (merged molecule, number of molecules merged)
    """
    merged_mol = Chem.RWMol()
    mols = rdk_read_molecule_files(inputs)
    index = 0
    for index, mol in enumerate(mols):
        if mol:
            merged_mol.InsertMol(mol)
        else:
            DmLog.emit_event("WARNING: could not process molecule", index)
    Chem.SanitizeMol(merged_mol)
    return merged_mol, index + 1


def rdk_mol_supplier(filename):
    """Generate a ForwardSDMolSupplier from a .sdf or .sdf.gz file."""
    if filename.endswith(".sdf"):
        return Chem.ForwardSDMolSupplier(filename)
    if filename.endswith(".sdf.gz"):
        return Chem.ForwardSDMolSupplier(gzip.open(filename, "rb"))
    raise ValueError(
        "Unsupported file type. Must be .sdf or .sdf.gz. Found " + filename
    )


def fragment(mol, mode):
    """Generate the largest fragment in the molecule, e.g. typically a
    desalting operation.

    :param mol: The molecule to fragment
    :param mode: The strategy for picking the largest fragment (mw or hac)
    :return: The largest fragment, with the original molecule's properties
    """
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol

    if mode == "hac":
        biggest_mol = max(frags, key=lambda frag: frag.GetNumHeavyAtoms())
    elif mode == "mw":
        biggest_mol = max(frags, key=Descriptors.MolWt)
    else:
        raise ValueError("Invalid fragment mode:", mode)

    # copy the properties across
    for name in mol.GetPropNames():
        biggest_mol.SetProp(name, mol.GetProp(name))

    # _Name is a magical property not included in GetPropNames
    if "_Name" in mol.GetPropNames():
        biggest_mol.SetProp("_Name", mol.GetProp("_Name"))

    return biggest_mol


def fragmentAndFingerprint(  # pylint: disable=invalid-name
    reader, mols, data, fps, descriptor, fragmentMethod="hac", outputFragment=False
):
    """Fragment the molecule if it has multiple fragments and generate
    fingerprints on the fragment.

    :param reader: MolSupplier from which to read the molecules
    :param mols: List to which the molecules are added
    :param data: List to which (mol_id, smiles, props) tuples are added
    :param fps: List to which the fingerprints are added
    :param descriptor: Function to generate the fingerprint from the molecule
    :param fragmentMethod: The fragmentation method for multi-fragment
        molecules (hac or mw)
    :param outputFragment: Whether to add the fragment or the original
        molecule to the mols list
    :return: The number of errors encountered
    """
    errors = 0
    count = 0
    while True:
        count += 1
        rec = reader.read()
        if not rec:
            break
        mol, _, mol_id, props = rec
        if not mol:
            log("Failed to read molecule", count)
            errors += 1
            continue

        frag = fragment(mol, fragmentMethod)
        fingerprint = descriptor(frag)
        if fingerprint:
            if outputFragment:
                mols.append(frag)
                data.append((mol_id, Chem.MolToSmiles(frag), props))
            else:
                mols.append(mol)
                data.append((mol_id, Chem.MolToSmiles(mol), props))
            fps.append(fingerprint)
    return errors


def check_molecules_are_3d(filename, num_to_check=10):
    """Check whether the first ``num_to_check`` molecules in a .sdf/.sdf.gz
    file have 3D conformers.
    """
    if filename.endswith(".gz"):
        suppl = Chem.ForwardSDMolSupplier(gzip.open(filename))
    else:
        suppl = Chem.ForwardSDMolSupplier(filename)

    count = 0
    for mol in suppl:
        conf = mol.GetConformer()
        if conf is None or not conf.Is3D():
            return False
        count += 1
        if count == num_to_check:
            break
    return True

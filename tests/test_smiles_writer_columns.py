import pathlib
import tempfile
import unittest

from rdkit import Chem

from rdkit_utils import create_reader, create_writer


def _write_one(**writer_kwargs):
    outfile = str(pathlib.Path(tempfile.mkdtemp()) / "out.smi")
    writer = create_writer(
        outfile, delimiter="\t", extra_field_names=[], calc_prop_names=["x"],
        **writer_kwargs,
    )
    writer.write("CCO", Chem.MolFromSmiles("CCO"), "mol1", [], [1.0])
    writer.close()
    return pathlib.Path(outfile).read_text(encoding="utf-8").strip()


class TestSmilesWriterMolColumn(unittest.TestCase):
    """mol_column=None must not survive into write(), which orders the ID and
    SMILES columns by comparing id_column with mol_column."""

    def test_none_mol_column_with_an_id_column_writes(self):
        # Previously raised TypeError: '<' not supported between 'int' and
        # 'NoneType'. Reachable from the CLI whenever a Job passes --id-column
        # and leaves --mol-column at its default.
        self.assertEqual(_write_one(id_column=1, mol_column=None), "CCO\tmol1\t1.0")

    def test_none_mol_column_without_an_id_column_writes(self):
        self.assertEqual(_write_one(id_column=None, mol_column=None), "CCO\t1.0")

    def test_id_column_zero_puts_the_id_first(self):
        self.assertEqual(_write_one(id_column=0, mol_column=None), "mol1\tCCO\t1.0")

    def test_explicit_mol_column_is_respected(self):
        self.assertEqual(_write_one(id_column=1, mol_column=0), "CCO\tmol1\t1.0")
        self.assertEqual(_write_one(id_column=0, mol_column=1), "mol1\tCCO\t1.0")

    def test_inference_matches_the_reader(self):
        # A writer built from the same arguments as its reader must agree with
        # it on where the molecule column sits.
        for id_column in (None, 0, 1, 2):
            with self.subTest(id_column=id_column):
                reader = create_reader(
                    str(pathlib.Path("test-data") / "mols.smi"),
                    delimiter="\t",
                    read_header=True,
                    id_column=id_column,
                    mol_column=None,
                )
                outfile = str(pathlib.Path(tempfile.mkdtemp()) / "out.smi")
                writer = create_writer(
                    outfile, delimiter="\t", extra_field_names=[],
                    calc_prop_names=[], id_column=id_column, mol_column=None,
                )
                self.assertEqual(writer.mol_column, reader.mol_column)
                writer.close()


if __name__ == "__main__":
    unittest.main()

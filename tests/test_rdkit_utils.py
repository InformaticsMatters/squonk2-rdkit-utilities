import pathlib
import tempfile
import unittest

from rdkit import Chem

from rdkit_utils import (
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

TEST_DATA = pathlib.Path("test-data")


class TestSmilesRoundTrip(unittest.TestCase):

    def test_read(self):
        reader = create_reader(
            str(TEST_DATA / "mols.smi"),
            delimiter="\t",
            read_header=True,
            id_column=1,
            mol_column=0,
        )
        records = []
        while True:
            rec = reader.read()
            if not rec:
                break
            records.append(rec)
        reader.close()

        self.assertEqual(len(records), 3)
        mol, smi, mol_id, props = records[0]
        self.assertIsNotNone(mol)
        self.assertEqual(mol_id, "mol1")
        self.assertEqual(Chem.CanonSmiles(smi), Chem.CanonSmiles("CCO"))
        self.assertEqual(props, ["ethanol"])

    def test_round_trip(self):
        reader = create_reader(
            str(TEST_DATA / "mols.smi"),
            delimiter="\t",
            read_header=True,
            id_column=1,
            mol_column=0,
        )
        extra_field_names = reader.get_extra_field_names()

        with tempfile.TemporaryDirectory() as tmp:
            outfile = str(pathlib.Path(tmp) / "out.smi")
            writer = create_writer(
                outfile,
                delimiter="\t",
                extra_field_names=extra_field_names,
                id_column=1,
                mol_column=0,
            )
            writer.write_header(["SMILES", "ID"] + extra_field_names)

            written_ids = []
            while True:
                rec = reader.read()
                if not rec:
                    break
                mol, smi, mol_id, props = rec
                writer.write(smi, mol, mol_id, props, [])
                written_ids.append(mol_id)
            reader.close()
            writer.close()

            reread = create_reader(
                outfile, delimiter="\t", read_header=True, id_column=1, mol_column=0
            )
            reread_ids = []
            while True:
                rec = reread.read()
                if not rec:
                    break
                reread_ids.append(rec[2])
            reread.close()

        self.assertEqual(written_ids, ["mol1", "mol2", "mol3"])
        self.assertEqual(reread_ids, written_ids)


class TestSdfRoundTrip(unittest.TestCase):

    def test_read(self):
        reader = create_reader(str(TEST_DATA / "mols.sdf"), id_column="_Name")
        mol, smi, mol_id, props = reader.read()
        reader.close()

        self.assertIsNotNone(mol)
        self.assertEqual(mol_id, "mol1")
        self.assertEqual(Chem.CanonSmiles(smi), Chem.CanonSmiles("CCO"))
        self.assertIn("ethanol", props)

    def test_round_trip(self):
        reader = create_reader(str(TEST_DATA / "mols.sdf"), id_column="_Name")
        extra_field_names = reader.get_extra_field_names()

        with tempfile.TemporaryDirectory() as tmp:
            outfile = str(pathlib.Path(tmp) / "out.sdf")
            writer = create_writer(
                outfile, calc_prop_names=["NEW_PROP"], id_column="_Name"
            )

            written_ids = []
            while True:
                rec = reader.read()
                if not rec:
                    break
                mol, smi, mol_id, props = rec
                writer.write(smi, mol, mol_id, props, ["value"])
                written_ids.append(mol_id)
            reader.close()
            writer.close()

            reread = create_reader(outfile, id_column="_Name")
            reread_ids = []
            reread_new_prop = []
            while True:
                rec = reread.read()
                if not rec:
                    break
                reread_ids.append(rec[2])
                reread_new_prop.append(rec[0].GetProp("NEW_PROP"))
            reread.close()

        self.assertEqual(written_ids, ["mol1", "mol2", "mol3"])
        self.assertEqual(reread_ids, written_ids)
        self.assertEqual(reread_new_prop, ["value"] * 3)


class TestGenerateHeaders(unittest.TestCase):

    def test_no_id_column(self):
        headers = generate_headers(0, None, None, ["field1"], ["calc1"], False)
        self.assertEqual(headers, ["SMILES", "field1", "calc1"])

    def test_sdf_id_column(self):
        headers = generate_headers(-1, "_Name", None, ["field1"], [], False)
        self.assertEqual(headers, ["SMILES", "ID", "field1"])

    def test_csv_id_column_with_omit_fields(self):
        headers = generate_headers(1, 1, None, ["SMILES", "ID"], ["calc1"], True)
        self.assertEqual(headers, ["SMILES", "ID", "calc1"])


class TestFragment(unittest.TestCase):

    def setUp(self):
        self.mol = Chem.MolFromSmiles("CCO.[Na+]")
        self.mol.SetProp("_Name", "sodium ethoxide salt")

    def test_fragment_hac(self):
        biggest = fragment(self.mol, "hac")
        self.assertEqual(Chem.CanonSmiles(Chem.MolToSmiles(biggest)), "CCO")
        self.assertEqual(biggest.GetProp("_Name"), "sodium ethoxide salt")

    def test_fragment_mw(self):
        biggest = fragment(self.mol, "mw")
        self.assertEqual(Chem.CanonSmiles(Chem.MolToSmiles(biggest)), "CCO")

    def test_single_fragment_returned_unchanged(self):
        mol = Chem.MolFromSmiles("CCO")
        self.assertIs(fragment(mol, "hac"), mol)

    def test_invalid_mode_raises(self):
        with self.assertRaises(ValueError):
            fragment(self.mol, "bogus")


class TestMolInspection(unittest.TestCase):

    def test_get_num_chiral_centers(self):
        # (R)/(S)-alanine has one chiral centre
        mol = Chem.MolFromSmiles("C[C@@H](N)C(=O)O")
        num_centers, undefined = get_num_chiral_centers(mol)
        self.assertEqual(num_centers, 1)
        self.assertEqual(undefined, 0)

    def test_get_num_sp3_centres(self):
        # ethanol: both carbons and the oxygen are sp3
        mol = Chem.MolFromSmiles("CCO")
        self.assertEqual(get_num_sp3_centres(mol), 3)

    def test_benzene_has_no_sp3_centres(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        self.assertEqual(get_num_sp3_centres(mol), 0)


class TestCheckMoleculesAre3d(unittest.TestCase):

    def test_2d_molecules_are_not_3d(self):
        self.assertFalse(check_molecules_are_3d(str(TEST_DATA / "mols.sdf")))

    def test_3d_molecules_are_3d(self):
        self.assertTrue(check_molecules_are_3d(str(TEST_DATA / "mols-3d.sdf")))


class TestMolFileHelpers(unittest.TestCase):

    def test_rdk_read_single_mol(self):
        mol = rdk_read_single_mol(str(TEST_DATA / "ethanol.mol"))
        self.assertEqual(Chem.CanonSmiles(Chem.MolToSmiles(mol)), "CCO")

    def test_rdk_read_mols(self):
        mols = rdk_read_mols(str(TEST_DATA / "ethanol.mol"))
        self.assertEqual(len(mols), 1)

    def test_rdk_read_molecule_files(self):
        mols = rdk_read_molecule_files(
            [f"{TEST_DATA / 'ethanol.mol'},{TEST_DATA / 'benzene.mol'}"]
        )
        self.assertEqual(len(mols), 2)

    def test_rdk_merge_mols(self):
        merged, count = rdk_merge_mols(
            [str(TEST_DATA / "ethanol.mol"), str(TEST_DATA / "benzene.mol")]
        )
        self.assertEqual(count, 2)
        self.assertEqual(merged.GetNumAtoms(), 3 + 6)


class TestFragmentAndFingerprint(unittest.TestCase):

    def test_fragments_and_fingerprints_each_molecule(self):
        reader = create_reader(str(TEST_DATA / "mols.sdf"), id_column="_Name")
        mols, data, fps = [], [], []
        errors = fragmentAndFingerprint(
            reader, mols, data, fps, lambda mol: mol.GetNumAtoms()
        )
        reader.close()

        self.assertEqual(errors, 0)
        self.assertEqual(len(mols), 3)
        self.assertEqual(len(fps), 3)
        self.assertEqual([mol_id for mol_id, _, _ in data], ["mol1", "mol2", "mol3"])


class TestSdfRecordGen(unittest.TestCase):

    def test_yields_one_record_per_molecule(self):
        with open(TEST_DATA / "mols.sdf", "rb") as handle:
            records = list(sdf_record_gen(handle))
        self.assertEqual(len(records), 3)
        for record in records:
            self.assertTrue(record.endswith("$$$$\n"))


class TestRdkMolSupplier(unittest.TestCase):

    def test_invalid_extension_raises(self):
        with self.assertRaises(ValueError):
            rdk_mol_supplier(str(TEST_DATA / "ethanol.mol"))


class TestUpdateChargeFlagInAtomBlock(unittest.TestCase):

    def test_matches_expected_output(self):
        molblock_in = (TEST_DATA / "molblock-a-in.sdf").read_text()
        expected = (TEST_DATA / "molblock-a-out.sdf").read_text()
        self.assertEqual(updateChargeFlagInAtomBlock(molblock_in), expected)


if __name__ == "__main__":
    unittest.main()

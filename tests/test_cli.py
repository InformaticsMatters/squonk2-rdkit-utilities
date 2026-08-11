import argparse
import pathlib
import unittest

from dm_job_utilities.utils import read_delimiter

from rdkit_utils import add_common_molecule_io_args, create_reader, str_or_int

TEST_DATA = pathlib.Path("test-data")


def _parser(**kwargs):
    parser = argparse.ArgumentParser()
    add_common_molecule_io_args(parser, **kwargs)
    return parser


class TestStrOrInt(unittest.TestCase):

    def test_int_like_value_returns_int(self):
        result = str_or_int("3")
        self.assertEqual(result, 3)
        self.assertIsInstance(result, int)

    def test_negative_int_like_value_returns_int(self):
        self.assertEqual(str_or_int("-1"), -1)

    def test_field_name_returns_original_string(self):
        self.assertEqual(str_or_int("SMILES"), "SMILES")

    def test_underscore_name_returns_original_string(self):
        self.assertEqual(str_or_int("_Name"), "_Name")


class TestAddCommonMoleculeIoArgs(unittest.TestCase):

    def test_parses_expected_namespace(self):
        args = _parser(include_y_column=True).parse_args(
            [
                "-i", "input.smi",
                "-o", "output.smi",
                "-d", "tab",
                "--id-column", "0",
                "--mol-column", "SMILES",
                "--y-column", "2",
                "--read-header",
                "--write-header",
                "--read-records", "50",
                "--omit-fields",
            ]
        )
        self.assertEqual(args.input, "input.smi")
        self.assertEqual(args.output, "output.smi")
        self.assertEqual(args.delimiter, "tab")
        self.assertEqual(args.id_column, 0)
        self.assertEqual(args.mol_column, "SMILES")
        self.assertEqual(args.y_column, 2)
        self.assertTrue(args.read_header)
        self.assertTrue(args.write_header)
        self.assertEqual(args.read_records, 50)
        self.assertTrue(args.omit_fields)

    def test_defaults_match_the_reader_signature(self):
        # mol_column and delimiter must default to None so that SmilesReader's
        # own inference applies; read_records matches create_reader().
        args = _parser().parse_args(["-i", "in.smi"])
        self.assertIsNone(args.output)
        self.assertIsNone(args.delimiter)
        self.assertIsNone(args.id_column)
        self.assertIsNone(args.mol_column)
        self.assertEqual(args.read_records, 100)
        self.assertFalse(args.read_header)
        self.assertFalse(args.write_header)
        self.assertFalse(args.omit_fields)

    def test_input_is_required(self):
        with self.assertRaises(SystemExit):
            _parser().parse_args([])

    def test_output_is_optional_by_default(self):
        args = _parser().parse_args(["-i", "in.smi"])
        self.assertIsNone(args.output)

    def test_output_can_be_required(self):
        with self.assertRaises(SystemExit):
            _parser(output_required=True).parse_args(["-i", "in.smi"])

    def test_output_default_is_applied(self):
        args = _parser(output_default="descriptors2d.smi").parse_args(
            ["-i", "in.smi"]
        )
        self.assertEqual(args.output, "descriptors2d.smi")

    def test_y_column_is_absent_unless_requested(self):
        args = _parser().parse_args(["-i", "in.smi"])
        self.assertFalse(hasattr(args, "y_column"))
        with self.assertRaises(SystemExit):
            _parser().parse_args(["-i", "in.smi", "--y-column", "2"])

    def test_y_column_defaults_to_none_when_included(self):
        args = _parser(include_y_column=True).parse_args(["-i", "in.smi"])
        self.assertIsNone(args.y_column)

    def test_returns_the_group_for_further_options(self):
        parser = argparse.ArgumentParser()
        group = add_common_molecule_io_args(parser)
        group.add_argument("--extra")
        args = parser.parse_args(["-i", "in.smi", "--extra", "value"])
        self.assertEqual(args.extra, "value")


class TestLegacyFlagAliases(unittest.TestCase):
    """--input/--output are retained so adopting the helper does not break
    existing Job manifests."""

    def test_input_is_an_alias_for_infile(self):
        args = _parser().parse_args(["--input", "in.smi"])
        self.assertEqual(args.input, "in.smi")

    def test_output_is_an_alias_for_outfile(self):
        args = _parser().parse_args(["-i", "in.smi", "--output", "out.smi"])
        self.assertEqual(args.output, "out.smi")

    def test_canonical_spellings_populate_the_same_destinations(self):
        args = _parser().parse_args(["--infile", "in.smi", "--outfile", "out.smi"])
        self.assertEqual(args.input, "in.smi")
        self.assertEqual(args.output, "out.smi")

    def test_short_forms_are_unchanged(self):
        args = _parser().parse_args(
            ["-i", "in.smi", "-o", "out.smi", "-k", "-d", "comma"]
        )
        self.assertEqual(args.input, "in.smi")
        self.assertEqual(args.output, "out.smi")
        self.assertTrue(args.omit_fields)
        self.assertEqual(args.delimiter, "comma")


class TestArgsFeedTheReaders(unittest.TestCase):
    """The point of the group: its namespace should drop straight into
    create_reader() without per-Job massaging."""

    def _read_all(self, argv):
        args = _parser().parse_args(argv)
        reader = create_reader(
            args.input,
            delimiter=read_delimiter(args.delimiter),
            read_header=args.read_header,
            id_column=args.id_column,
            mol_column=args.mol_column,
            read_records=args.read_records,
        )
        records = []
        while True:
            record = reader.read()
            if record is None:
                break
            records.append(record)
        return reader, records

    def test_defaults_drive_create_reader(self):
        reader, records = self._read_all(
            ["-i", str(TEST_DATA / "mols.smi"), "-d", "tab", "--read-header"]
        )
        # mol_column was left at its None default, so SmilesReader inferred it
        self.assertEqual(reader.mol_column, 0)
        self.assertIsNone(reader.id_column)
        self.assertEqual(reader.field_names, ["SMILES", "ID", "NAME"])
        self.assertEqual(len(records), 3)

    def test_explicit_id_column_shifts_the_inferred_mol_column(self):
        reader, records = self._read_all(
            [
                "-i", str(TEST_DATA / "mols.smi"),
                "-d", "tab",
                "--read-header",
                "--id-column", "0",
            ]
        )
        self.assertEqual(reader.id_column, 0)
        self.assertEqual(reader.mol_column, 1)
        self.assertEqual(len(records), 3)


if __name__ == "__main__":
    unittest.main()

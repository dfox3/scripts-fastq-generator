import unittest
import pytest
import scripts_fastq_generator.fastq_generator as fg


def test_generate_random_fastq_set():
    assert True


def test_generate_record():
    assert True


class TestGenerateTitle(unittest.TestCase):
    def test_title_clean(self):
        desc = "READ67"
        self.assertEqual(fg.generate_title(desc), "READ67")



def test_generate_q_string():
    assert True


def test__generate_read_from_sequence():
    assert True


def test__rand_read_from_seq():
    assert True


def test_get_sequences_from_bed():
    assert True


def test_add_error():
    assert True


def test__rand_error_base():
    assert True


def test_load_sequences_from_fasta():
    assert True


class TestReverseComplement(unittest.TestCase):
    def test_reverse_complement_clean(self):
        sequence = "AATTCCGG"
        self.assertEqual(fg.reverse_complement(sequence), "CCGGAATT")

    def test_reverse_complement_palindrome(self):
        sequence = "ATCGAT"
        self.assertEqual(fg.reverse_complement(sequence), "ATCGAT")

    def test_reverse_complement_wrong_types(self):
        sequence = 234345234.09385
        with pytest.raises(TypeError) as pytest_wrapped_e:
            fg.reverse_complement(sequence)
        self.assertEqual(pytest_wrapped_e.type, TypeError)

        sequence = None
        with pytest.raises(TypeError) as pytest_wrapped_e:
            fg.reverse_complement(sequence)
        self.assertEqual(pytest_wrapped_e.type, TypeError)

        sequence = ["A", "C", "T", "G"]
        with pytest.raises(AttributeError) as pytest_wrapped_e:
            fg.reverse_complement(sequence)
        self.assertEqual(pytest_wrapped_e.type, AttributeError)

    def test_reverse_complement_nonsense_string(self):
        sequence = "ihd939q84a;sda"
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            fg.reverse_complement(sequence)
        self.assertEqual(pytest_wrapped_e.type, SystemExit)
        self.assertEqual(pytest_wrapped_e.value.code, 2)


def test_write_gzip_string():
    assert True

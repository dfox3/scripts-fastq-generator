import unittest
import pytest
import scripts_fastq_generator.fastq_generator as fg
import subprocess
import sys


class TestMain(unittest.TestCase):
    def test_all_t(self):
        result = subprocess.run(
            [
                "python",
                "scripts_fastq_generator/fastq_generator.py",
                "-o",
                "test_out/dumpster",
                "-i",
                "tests/fixtures/all_T.fasta",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0)

    def test_all_wrong_ref(self):
        result = subprocess.run(
            [
                "python",
                "scripts_fastq_generator/fastq_generator.py",
                "-o",
                "test_out/dumpster",
                "-b",
                "tests/fixtures/all_T.fasta",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)

    def test_all_no_output(self):
        result = subprocess.run(
            [
                "python",
                "scripts_fastq_generator/fastq_generator.py",
                "-i",
                "tests/fixtures/all_T.fasta",
            ],
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)

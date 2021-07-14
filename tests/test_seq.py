import pytest
import unittest
import alhambra.seq as seq
import random


class TestSeq(unittest.TestCase):
    seq1 = "aggagtacc"
    seq2 = "aggn"
    seq3 = "gngc"
    seq4 = "hggc"
    defbases = frozenset(("a", "g", "t", "c", "A", "G", "C", "T"))
    ambbases = frozenset(seq._L_TO_N.keys()) - defbases
    allbases = frozenset(seq._L_TO_N.keys())

    def test_merge_correct(self):
        # same sequence
        self.assertEqual(seq.merge(self.seq1, self.seq1), self.seq1)
        # constraints on both sides
        self.assertEqual(seq.merge(self.seq2, self.seq4), "aggc")
        self.assertEqual(seq.merge(self.seq4, self.seq2), "aggc")

    def test_merge_conflicts(self):
        with self.assertRaises(seq.MergeConflictError):
            seq.merge(self.seq1, self.seq2)
        with self.assertRaises(seq.MergeConflictError):
            seq.merge(self.seq2, self.seq3)
        with self.assertRaises(seq.MergeConflictError):
            seq.merge(self.seq3, self.seq4)

    def test_is_null(self):
        self.assertTrue(seq.is_null(""))
        self.assertTrue(seq.is_null("nnnnnnnnnnnnnnnnNNNNNNN"))
        self.assertFalse(seq.is_null("nnagaagannna"))
        self.assertFalse(seq.is_null("nnnnnnnynnnn"))

    def test_is_definite(self):
        for _ in range(0, 100):
            self.assertTrue(
                seq.is_definite("".join(random.choices(tuple(self.defbases), k=50)))
            )
            self.assertFalse(
                seq.is_definite("".join(random.choices(tuple(self.ambbases), k=50)))
            )
            self.assertFalse(
                seq.is_definite(
                    "".join(
                        random.choices(
                            tuple(frozenset.union(self.defbases, self.ambbases)), k=50
                        )
                    )
                    + "b"
                )
            )
            # FIXME: this last test needs that extra b to ensure that random doesn't choose only
            # definite bases.  Test coverage is hindered.

    def test_space(self):
        self.assertTrue(seq.is_null("Nn nn  nn\tnn\nnn"))
        self.assertTrue(seq.is_null("       \t\n  \n"))
        self.assertFalse(seq.is_definite("    \t\n  \n"))
        self.assertFalse(seq.is_null("nnnnnnnn      \t a \n nnnnnn"))
        self.assertFalse(seq.is_definite("ag  n n\t\nn"))
        self.assertTrue(seq.is_definite("ag  c t\t\ntt"))

    def test_length(self):
        n = random.randint(0, 100)
        s = list(
            "".join(random.choices(tuple(self.allbases), k=n))
            + "".join(random.choices(tuple(seq._PUNC), k=random.randint(0, 100)))
        )
        random.shuffle(s)
        self.assertEqual(seq.dna_length("".join(s)), n)
        self.assertEqual(seq.dna_length("\n\t"), 0)
        self.assertEqual(seq.dna_length(""), 0)

    def test_count_ambiguous(self):
        for _ in range(0, 100):
            n = random.randint(0, 1000)
            s = list(
                "".join(random.choices(tuple(self.ambbases), k=n))
                + "".join(
                    random.choices(tuple(self.defbases), k=random.randint(0, 1000))
                )
                + "".join(random.choices(tuple(seq._PUNC), k=random.randint(0, 1000)))
            )
            random.shuffle(s)
            self.assertEqual(seq.count_ambiguous("".join(s)), n)
        # FIXME: deal with uppercase.

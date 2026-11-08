import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location('counts', Path(__file__).with_name('04_consensus_counts.py'))
counts = importlib.util.module_from_spec(spec); spec.loader.exec_module(counts)


class CountChecks(unittest.TestCase):
    def test_overlap_counts_against_exhaustive_pairs(self):
        rng = np.random.default_rng(20260908)
        starts = rng.integers(0, 10000, size=1000)
        ends = starts+rng.integers(1, 2000, size=1000)
        region_starts = np.arange(0, 10000, 501)
        region_ends = region_starts+501
        expected = np.array([sum((starts < b) & (ends > a)) for a,b in zip(region_starts,region_ends)])
        actual = counts.interval_counts(starts,ends,region_starts,region_ends)
        np.testing.assert_array_equal(actual,expected)

    def test_frip_counts_long_fragments_only_once(self):
        starts = np.array([0,100,200,50,400])
        ends = np.array([300,150,220,100,450])
        regions_start = np.array([100,200]); regions_end = np.array([150,250])
        per_peak = counts.interval_counts(starts,ends,regions_start,regions_end)
        self.assertEqual(per_peak.sum(),4)
        self.assertEqual(counts.fragments_in_union(starts,ends,regions_start,regions_end),3)


if __name__ == '__main__': unittest.main(verbosity=2)

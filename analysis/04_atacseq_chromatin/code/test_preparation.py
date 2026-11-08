"""Small biological edge-case checks before scanning the real BAMs."""
from collections import Counter
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pysam
import pyBigWig

spec = importlib.util.spec_from_file_location('prepare', Path(__file__).with_name('02_prepare.py'))
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class PreparationChecks(unittest.TestCase):
    def pair(self, header, name, first, second, width=50, mapq=60, duplicate=False, chrom='chr1'):
        records = []
        span = second + width - first
        for index, position in enumerate([first, second]):
            read = pysam.AlignedSegment(header)
            read.query_name = name; read.query_sequence = 'A' * width
            read.flag = (99 if index == 0 else 147) | (1024 if duplicate else 0)
            read.reference_id = header.get_tid(chrom); read.reference_start = position
            read.mapping_quality = 60 if index == 0 else mapq
            read.cigarstring = f'{width}M'
            read.next_reference_id = read.reference_id
            read.next_reference_start = second if index == 0 else first
            read.template_length = span if index == 0 else -span
            records.append(read)
        return records

    def test_whole_pair_filtering(self):
        header = pysam.AlignmentHeader.from_dict({'HD': {'SO': 'coordinate'},
          'SQ': [{'SN': 'chr1', 'LN': 1000}, {'SN': 'chrM', 'LN': 1000}]})
        records = []
        records += self.pair(header, 'goodA', 100, 150)
        records += self.pair(header, 'goodB', 110, 160)
        records += self.pair(header, 'blacklisted_gap', 210, 280, width=20)
        records += self.pair(header, 'one_low_mapq', 400, 450, mapq=10)
        records += self.pair(header, 'orphan', 600, 650)[:1]
        records += self.pair(header, 'duplicate', 700, 750, duplicate=True)
        records += self.pair(header, 'mitochondrial', 100, 150, chrom='chrM')
        with tempfile.TemporaryDirectory() as tmp:
            bam_path = Path(tmp)/'test.bam'
            with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as bam:
                for read in sorted(records, key=lambda r: (r.reference_id, r.reference_start)):
                    bam.write(read)
            metrics = Counter()
            with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
                pairs = list(prepare.eligible_pairs(bam, {'chr1'}, {'chr1': ([250], [260])}, metrics))
        self.assertEqual({p[0].query_name for p in pairs}, {'goodA', 'goodB'})
        self.assertEqual(metrics['accepted_fragments'], 2)
        self.assertEqual(metrics['low_mapq_pairs'], 1)
        self.assertEqual(metrics['blacklist_pairs'], 1)
        self.assertEqual(metrics['duplicate_pairs'], 1)
        self.assertEqual(metrics['unmatched_primary_records'], 1)
        self.assertEqual(metrics['noncanonical_records'], 2)

    def test_half_open_blacklist_boundaries(self):
        blacklist = {'chr1': ([100, 200], [150, 250])}
        self.assertFalse(prepare.overlaps(blacklist, 'chr1', 50, 100))
        self.assertFalse(prepare.overlaps(blacklist, 'chr1', 150, 200))
        self.assertTrue(prepare.overlaps(blacklist, 'chr1', 99, 101))
        self.assertTrue(prepare.overlaps(blacklist, 'chr1', 50, 300))

    def test_track_preserves_partial_bin_coverage(self):
        fragments = {'chr1': (np.array([10,20], dtype=np.uint32), np.array([40,30], dtype=np.uint32))}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'test.bw'
            prepare.write_coverage(path, fragments, {'chr1': 60}, denominator=1e6)
            with pyBigWig.open(str(path)) as bw:
                self.assertAlmostEqual(bw.values('chr1', 0, 1)[0], 0.8, places=6)
                self.assertAlmostEqual(bw.values('chr1', 25, 26)[0], 0.8, places=6)


if __name__ == '__main__':
    unittest.main(verbosity=2)

import pysam
import pytest
from collections import Counter
from filter_junction_chain_sr_or_lr_support_26_08_26 import (count_junctions, junction_support, parse_gtf_exons,
                                                             build_junctions, min_junction_coverage)

# ---------- build_junctions: the off-by-one ----------


def test_junction_coords():
    exons = {"t1": [(1000, 1100), (1601, 1700)]}
    chroms = {"t1": "chr1"}
    txj, allj = build_junctions(exons, chroms)
    assert txj["t1"] == [("chr1", 1100, 1600)]
    assert allj == {("chr1", 1100, 1600)}


def test_multi_intron_and_shared():
    exons = {"t1": [(100, 200), (301, 400), (501, 600)],
             "t2": [(100, 200), (301, 400)]}
    chroms = {"t1": "chr1", "t2": "chr1"}
    txj, allj = build_junctions(exons, chroms)
    assert len(txj["t1"]) == 2
    assert len(allj) == 2          # shared junction deduplicated


def test_mono_exon_has_no_junctions():
    txj, allj = build_junctions({"t1": [(100, 200)]}, {"t1": "chr1"})
    assert txj["t1"] == [] and allj == set()

# ---------- parse_gtf_exons ----------


def _gtf(tmp_path, lines):
    p = tmp_path / "x.gtf"
    p.write_text("".join(lines))
    return str(p)


def test_parse_sorts_and_skips(tmp_path):
    lines = ["#comment\n",
             'chr1\ts\ttranscript\t1000\t1700\t.\t-\t.\ttranscript_id "t1";\n',
             'chr1\ts\texon\t1601\t1700\t.\t-\t.\ttranscript_id "t1";\n',
             'chr1\ts\texon\t1000\t1100\t.\t-\t.\ttranscript_id "t1";\n']
    exons, chroms = parse_gtf_exons(_gtf(tmp_path, lines))
    assert exons["t1"] == [(1000, 1100), (1601, 1700)
                           ]   # minus-strand reordered
    assert chroms["t1"] == "chr1"


def test_overlapping_exons_raise(tmp_path):
    lines = ['chr1\ts\texon\t1000\t1100\t.\t+\t.\ttranscript_id "t1";\n',
             'chr1\ts\texon\t1050\t1200\t.\t+\t.\ttranscript_id "t1";\n']
    with pytest.raises(ValueError):
        parse_gtf_exons(_gtf(tmp_path, lines))

# ---------- junction_support: the tolerance ----------


def test_rescue_and_limits():
    allj = {("20", 1000, 2000)}
    sup = junction_support(Counter({("20", 1001, 2000): 3}), allj, tol=2)
    assert sup[("20", 1000, 2000)] == (0, 3)          # rescued

    sup = junction_support(Counter({("20", 1003, 2000): 3}), allj, tol=2)
    assert sup[("20", 1000, 2000)] == (0, 0)          # 3bp away, not rescued


def test_within_never_below_exact():
    allj = {("20", 1000, 2000)}
    c = Counter({("20", 1000, 2000): 50, ("20", 1001, 2000): 7})
    ex, wi = junction_support(c, allj, tol=2)[("20", 1000, 2000)]
    assert ex == 50 and wi == 57 and wi >= ex

# ---------- min_junction_coverage ----------


def test_min_and_mono():
    txj = {"t1": [("c", 1, 2), ("c", 3, 4)], "t2": []}
    sup = {("c", 1, 2): (10, 12), ("c", 3, 4): (0, 3)}
    mc = min_junction_coverage(txj, sup)
    assert mc["t1"] == (0, 3)         # weakest junction wins
    assert mc["t2"] == (None, None)   # mono-exonic

# ---------- count_junctions: the CIGAR walk ----------


@pytest.fixture
def bam(tmp_path):
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 100000, "SN": "chr1"}]}
    reads = [("plain",        1000, "100M500N100M",            60, 0),
             ("del_near_junc", 1000, "10S100M500N80M2D20M5S",   60, 0),
             ("insertion",    4000, "50M3I50M200N60M",         60, 0),
             ("tiny_anchor",  2000, "100M300N4M",              60, 0),
             ("short_intron", 3000, "50M20N50M",               60, 0),
             ("unspliced",    5000, "150M",                    60, 0),
             ("low_mapq",     7000, "80M400N80M",               0, 0),
             ("secondary",    1000, "100M500N100M",            60, 256)]
    p = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(p, "wb", header=hdr) as out:
        for name, start, cig, mapq, flag in sorted(reads, key=lambda r: r[1]):
            a = pysam.AlignedSegment()
            a.query_name, a.reference_id, a.reference_start = name, 0, start
            a.mapping_quality, a.cigarstring, a.flag = mapq, cig, flag
            n = sum(l for op, l in a.cigartuples if op in {0, 1, 4, 7, 8})
            a.query_sequence = "A" * n
            out.write(a)
    return p


def test_cigar_walk(bam):
    c = count_junctions(bam)
    # plain + del_near_junc, D kept coords right
    assert c[("chr1", 1100, 1600)] == 2
    # insertion did not shift ref pointer
    assert c[("chr1", 4100, 4300)] == 1
    assert ("chr1", 2100, 2400) not in c     # 4bp anchor rejected
    assert ("chr1", 3050, 3070) not in c     # 20bp "intron" rejected
    assert ("chr1", 7100, 7500) not in c     # mapq 0 rejected
    # unspliced + secondary contribute nothing
    assert len(c) == 2


def test_end_to_end(bam, tmp_path):
    lines = ['chr1\ts\texon\t1000\t1100\t.\t+\t.\ttranscript_id "t1";\n',
             'chr1\ts\texon\t1601\t1700\t.\t+\t.\ttranscript_id "t1";\n',
             'chr1\ts\texon\t9000\t9100\t.\t+\t.\ttranscript_id "t2";\n',
             'chr1\ts\texon\t9601\t9700\t.\t+\t.\ttranscript_id "t2";\n']
    exons, chroms = parse_gtf_exons(_gtf(tmp_path, lines))
    txj, allj = build_junctions(exons, chroms)
    mc = min_junction_coverage(txj, junction_support(
        count_junctions(bam), allj, tol=2))
    assert mc["t1"] == (2, 2)        # supported
    assert mc["t2"] == (0, 0)        # not supported

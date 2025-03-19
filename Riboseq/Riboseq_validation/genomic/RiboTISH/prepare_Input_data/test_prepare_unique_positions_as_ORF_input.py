"""
Script Name: test_prepare_unique_positions_as_ORF_input.py
Description: ...
.
Usage: python test_prepare_unique_positions_as_ORF_input.py
run tests on the functionality of the Region calculation functions for RiboTISH
in transcriptomic coordinates
"""

from prepare_unique_positions_as_ORF_input import create_region_for_RiboTISH

# test the supplied function
assert create_region_for_RiboTISH(0,  4,  9) == [3, 9]
assert create_region_for_RiboTISH(0,  4,  10) == [3, 12]
assert create_region_for_RiboTISH(1,  3,  10) == [1, 10]

assert create_region_for_RiboTISH(1,  3,  10) == [1, 10]
assert create_region_for_RiboTISH(1,  4,  10) == [4, 10]
assert create_region_for_RiboTISH(1,  3,  11) == [1, 13]

assert create_region_for_RiboTISH(2,  5,  11) == [5, 11]
assert create_region_for_RiboTISH(2,  5,  12) == [5, 14]
assert create_region_for_RiboTISH(2,  4,  11) == [2, 11]
assert create_region_for_RiboTISH(2,  4,  12) == [2, 14]
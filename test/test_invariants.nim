import std/unittest

import ../src/links
import ../src/invariants

test "linking_matrix":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.linking_matrix() == @[@[0]]
    let l2a1_0 = link_from_PD_code(@[[4, 1, 3, 2], [2, 3, 1, 4]])
    check l2a1_0.linking_matrix() == @[@[0, -1], @[-1, 0]]
    let l2a1_1 = link_from_PD_code(@[[4, 2, 3, 1], [2, 4, 1, 3]])
    check l2a1_1.linking_matrix() == @[@[0, 1], @[1, 0]]

test "colorability_matrix":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.colorability_matrix() == @[@[-1, -1, 2], @[2, -1, -1], @[-1, 2, -1]]
    let hopf_positive = link_from_PD_code(@[[2, 1, 3, 0], [1, 2, 0, 3]])
    check hopf_positive.colorability_matrix() == @[@[-2, 2], @[2, -2]]
    let hopf_negative = link_from_PD_code(@[[0, 2, 1, 3], [3, 1, 2, 0]])
    check hopf_negative.colorability_matrix() == @[@[-2, 2], @[2, -2]]

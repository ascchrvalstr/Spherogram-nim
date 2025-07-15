import std/unittest
import std/sets

import ../src/links
import ../src/seifert

test "seifert_circles":
    var fig8 = link_from_PD_code(@[[1, 7, 2, 6], [5, 3, 6, 2], [7, 4, 0, 5], [3, 0, 4, 1]])
    check fig8.seifert_circles() == [@[(0, 1), (2, 3), (1, 1), (3, 3)], @[(0, 2), (1, 2)], @[(2, 2), (3, 2)]]

test "admissible_moves":
    var fig8 = link_from_PD_code(@[[1, 7, 2, 6], [5, 3, 6, 2], [7, 4, 0, 5], [3, 0, 4, 1]])
    check fig8.admissible_moves() == (@[((0, 0), (3, 2)), ((3, 2), (2, 0)), ((0, 1), (1, 3)), ((2, 3), (1, 3)), ((0, 3), (1, 1)), ((0, 3), (3, 3)), ((1, 0), (2, 2)), ((2, 2), (3, 0))], @[(0, 2), (2, 0), (0, 1), (0, 1), (1, 0), (1, 0), (0, 2), (2, 0)])
    # check cs_to_seifert_circle == @[[0, 0, 1, 1], [0, 0, 1, 1], [0, 2, 2, 0], [0, 2, 2, 0]]

test "seifert_tree":
    var K5a2 = link_from_PD_code(@[[7, 3, 8, 2], [9, 5, 0, 4], [1, 7, 2, 6], [3, 9, 4, 8], [5, 1, 6, 0]])
    check K5a2.seifert_tree() == @[[toHashSet([0, 1]), toHashSet([0])], [toHashSet([1]), toHashSet([0, 1])]]
    # check circles == @[@[(0, 1), (3, 1), (1, 1), (4, 1), (2, 1)], @[(0, 2), (3, 2), (1, 2), (4, 2), (2, 2)]]
    # check strand_to_circle == @[[-1, 0, 1, -1], [-1, 0, 1, -1], [-1, 0, 1, -1], [-1, 0, 1, -1], [-1, 0, 1, -1]]
    # check (under_circle, over_circle) == (0, 1)

    var fig8 = link_from_PD_code(@[[7, 4, 0, 5], [3, 0, 4, 1], [1, 7, 2, 6], [5, 3, 6, 2]])
    check fig8.seifert_tree() == @[[toHashSet([0, 1]), toHashSet([0])],
                                   [toHashSet([1, 2]), toHashSet([0, 1])],
                                   [toHashSet([2]), toHashSet([1, 2])]]

test "remove_admissible_move":
    var knot52 = link_from_PD_code(@[[4, 0, 5, 9], [0, 6, 1, 5], [8, 2, 9, 1], [2, 8, 3, 7], [6, 4, 7, 3]])
    check knot52.remove_admissible_move()
    # check moves == @[((0, 1), (1, 3)), ((0, 3), (2, 1)), ((0, 3), (3, 3)), ((0, 3), (4, 1)), ((2, 1), (3, 3)), ((2, 1), (4, 1)), ((3, 3), (4, 1)), ((1, 1), (4, 3)), ((1, 1), (3, 1)), ((1, 1), (2, 3)), ((4, 3), (3, 1)), ((4, 3), (2, 3)), ((3, 1), (2, 3))]
    # check circle_pairs == @[(0, 1), (1, 2), (1, 3), (1, 0), (2, 3), (2, 0), (3, 0), (0, 3), (0, 2), (0, 1), (3, 2), (3, 1), (2, 1)]
    # check tree == @[[{3, 0, 2, 1}, {0}], [{1}, {3, 0, 2, 1}], [{3, 0, 2, 1}, {2}], [{3}, {3, 0, 2, 1}]]
    # check (e1_index, e2_index) == (0, 2)
    # check n == 5
    # check pair == (2, 0)
    # check cs1 == (2, 1)
    # check cs2 == (4, 1)

test "isotope_to_braid and is_chain":
    var knot52 = link_from_PD_code(@[[4, 0, 5, 9], [0, 6, 1, 5], [8, 2, 9, 1], [2, 8, 3, 7], [6, 4, 7, 3]])
    isotope_to_braid(knot52)
    check is_chain(seifert_tree(knot52))

    var knot61 = link_from_PD_code(@[[6, 11, 7, 0], [0, 5, 1, 6], [10, 2, 11, 1], [2, 10, 3, 9], [8, 4, 9, 3], [4, 8, 5, 7]])
    isotope_to_braid(knot61)
    check is_chain(seifert_tree(knot61))

test "straighten_arrows":
    var arrows = @[[0, 0, 0, 0], [1, 1, 0, 0], [2, 2, 0, 0]]
    straighten_arrows(arrows)
    check arrows == @[[0, 0, 0, 0], [1, 1, 0, 0], [2, 2, 0, 0]]

    arrows = @[[0, 0, 0, 1], [1, 2, 0, 1], [1, 0, 1, 0], [3, 1, 1, 0]]
    straighten_arrows(arrows)
    check arrows == @[[0, 0, 0, 1], [2, 2, 0, 1], [1, 1, 1, 0], [3, 3, 1, 0]]

    arrows = @[[0, 0, 0, 0], [1, 4, 0, 1], [1, 0, 1, 1], [2, 2, 1, 1], [3, 3, 1, 1], [5, 4, 1, 1], [6, 6, 1, 1], [1, 0, 2, 0], [5, 1, 2, 1]]
    straighten_arrows(arrows)
    check arrows == @[[0, 0, 0, 0], [5, 5, 0, 1], [1, 1, 1, 1], [3, 3, 1, 1], [4, 4, 1, 1], [6, 6, 1, 1], [8, 8, 1, 1], [2, 2, 2, 0], [7, 7, 2, 1]]

test "braid_arrows":
    var trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check braid_arrows(trefoil) == @[[0, 0, 0], [1, 0, 0], [2, 0, 0]]
    # check circles == @[@[(0, 2), (1, 2), (2, 2)], @[(0, 3), (1, 3), (2, 3)]]
    # check tree == @[[{0, 1}, {0}], [{1}, {0, 1}]]
    # check tails == @[{0, 1}, {1}]
    # check heads == {{0}, {0, 1}}
    # check start == 1
    # check ordered_strands == @[@[(0, 3), (1, 3), (2, 3)], @[(0, 2), (1, 2), (2, 2)]]
    # check ordered_strands == @[@[(0, 3), (1, 3), (2, 3)], @[(0, 2), (1, 2), (2, 2)]]
    # check positions_in_next_strand == @[{0: (0, 0), 1: (1, 0), 2: (2, 0)}]
    # check arrows == @[[2, 2, 0, 0], [0, 0, 0, 0], [1, 1, 0, 0]]
    # check arrows == @[[0, 0, 0, 0], [1, 1, 0, 0], [2, 2, 0, 0]]

test "braid_word":
    var trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check braid_word(trefoil) == [-1, -1, -1]

    var fig8 = link_from_PD_code(@[[1, 7, 2, 6], [5, 3, 6, 2], [7, 4, 0, 5], [3, 0, 4, 1]])
    check braid_word(fig8) == [1, -2, 1, -2]

    var knot51 = link_from_PD_code(@[[9, 4, 0, 5], [5, 0, 6, 1], [1, 6, 2, 7], [7, 2, 8, 3], [3, 8, 4, 9]])
    check braid_word(knot51) == [-1, -1, -1, -1, -1]

    var knot52 = link_from_PD_code(@[[4, 0, 5, 9], [0, 6, 1, 5], [8, 2, 9, 1], [2, 8, 3, 7], [6, 4, 7, 3]])
    check braid_word(knot52) == [-1, 2, 3, 2, 2, 1, 2, -3, 2]

test "seifert_matrix":
    var trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.seifert_matrix()[0] == @[@[-1, 0], @[1, -1]]

    var hopf_positive = link_from_PD_code(@[[2, 1, 3, 0], [1, 2, 0, 3]])
    check hopf_positive.seifert_matrix()[0] == @[@[1]]

    var hopf_negative = link_from_PD_code(@[[0, 2, 1, 3], [3, 1, 2, 0]])
    check hopf_negative.seifert_matrix()[0] == @[@[-1]]

    var fig8 = link_from_PD_code(@[[7, 4, 0, 5], [3, 0, 4, 1], [1, 7, 2, 6], [5, 3, 6, 2]])
    check fig8.seifert_matrix()[0] == @[@[-1, 0], @[1, 1]]

    var knot51 = link_from_PD_code(@[[9, 4, 0, 5], [5, 0, 6, 1], [1, 6, 2, 7], [7, 2, 8, 3], [3, 8, 4, 9]])
    check knot51.seifert_matrix()[0] == @[@[-1, 0, 0, 0], @[1, -1, 0, 0], @[0, 1, -1, 0], @[0, 0, 1, -1]]

    var knot52 = link_from_PD_code(@[[4, 0, 5, 9], [0, 6, 1, 5], [8, 2, 9, 1], [2, 8, 3, 7], [6, 4, 7, 3]])
    check knot52.seifert_matrix()[0] == @[@[0, 0, 0, 0, 0, 0], @[1, 1, -1, 0, 0, 0], @[0, 0, 1, -1, 0, 0], @[0, 0, 0, 1, -1, 0], @[-1, 0, 0, 0, 1, 0], @[0, 0, 0, 1, 0, 0]]

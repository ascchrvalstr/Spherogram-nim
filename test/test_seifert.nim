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

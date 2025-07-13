import std/unittest

import ../src/links
import ../src/twist

test "make_twist_regions_consistent":
    var empty_link = link_from_PD_code(@[])
    check make_twist_regions_consistent(empty_link) == 0
    check empty_link.crossings == []
    check empty_link.signs == []
    check empty_link.link_components == []

    var one_vertex_unknot = link_from_PD_code(@[[0, 0, 1, 1]])
    check make_twist_regions_consistent(one_vertex_unknot) == 0
    check one_vertex_unknot.crossings == [[1, 0, 3, 2]]
    check one_vertex_unknot.signs == [1]
    check one_vertex_unknot.link_components == [newLinkComponent[int](0, 1)]

    var unlinked_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 3, 1, 2]], @[1, -1])
    check make_twist_regions_consistent(unlinked_hopf) == 1
    check unlinked_hopf.crossings == [[5, 4, 7, 6], [1, 0, 3, 2]]
    check unlinked_hopf.signs == [1, 1]
    check unlinked_hopf.link_components == [newLinkComponent[int](0, 1), newLinkComponent[int](0, 2)]

    var two_crossing_switches = link_from_PD_code(@[[5, 2, 0, 3], [4, 4, 5, 3], [0, 2, 1, 1]])
    check make_twist_regions_consistent(two_crossing_switches) == 2
    check two_crossing_switches.crossings == [[7, 10, 9, 4], [3, 6, 5, 0], [11, 2, 1, 8]]
    check two_crossing_switches.signs == [-1, -1, -1]
    check two_crossing_switches.link_components == [newLinkComponent[int](0, 2)]

    var one_leads_to_two = link_from_PD_code(@[[5, 5, 0, 4], [0, 3, 1, 4], [2, 1, 3, 2]])
    check make_twist_regions_consistent(one_leads_to_two) == 2
    check one_leads_to_two.crossings == [[1, 0, 7, 6], [9, 8, 3, 2], [5, 4, 11, 10]]
    check one_leads_to_two.signs == [1, 1, 1]
    check one_leads_to_two.link_components == [newLinkComponent[int](0, 1)]

    var one_kills_off_the_other = link_from_PD_code(@[[5, 5, 0, 4], [0, 3, 1, 4], [1, 3, 2, 2]])
    check make_twist_regions_consistent(one_kills_off_the_other) == 1
    check one_kills_off_the_other.crossings == [[1, 0, 7, 6], [9, 8, 3, 2], [5, 4, 11, 10]]
    check one_kills_off_the_other.signs == [1, 1, 1]
    check one_kills_off_the_other.link_components == [newLinkComponent[int](0, 1)]

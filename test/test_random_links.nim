import std/unittest

import ../src/links
import ../src/random_links

test "random_link_internal":
    for num_crossings in 1 .. 10:
        var link = random_link_internal[int](num_crossings)
        check link.crossings.len == num_crossings

test "longest_components":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check longest_components(trefoil, 0) == []
    check longest_components(trefoil, 1) == [0]
    # check components_by_self_crossings == @[(3, 0)]

    let unlinked_hopf = link_from_PD_code(@[[0, 4, 1, 5], [3, 4, 0, 5], [1, 3, 2, 2]], @[-1, 1, 1])
    check longest_components(unlinked_hopf, 0) == []
    check longest_components(unlinked_hopf, 1) == [0]
    # check components_by_self_crossings == @[(1, 0), (0, 1)]
    check longest_components(unlinked_hopf, 2) == [0, 1]

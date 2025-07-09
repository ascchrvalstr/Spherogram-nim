import std/unittest
import std/options

import ../src/links
import ../src/random_links
import ../src/simplify

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

test "random_link":
    for num_crossings in 10 .. 15:
        var links = random_link(num_crossings, simplification_mode = none(SimplificationMode), prime_decomposition = false)
        check links.len == 1
        var link = links[0]
        check link.crossings.len == num_crossings
        links = random_link(num_crossings, simplification_mode = none(SimplificationMode))
        check links.len == 1
        link = links[0]
        check link.crossings.len <= num_crossings
        links = random_link(num_crossings, prime_decomposition = false)
        check links.len == 1
        link = links[0]
        check link.crossings.len <= num_crossings
        links = random_link(num_crossings, some(3), true, simplification_mode = none(SimplificationMode), prime_decomposition = false)
        check links.len == 1
        link = links[0]
        check link.crossings.len == num_crossings
        check link.link_components.len == 3
        links = random_link(num_crossings, some(3), simplification_mode = none(SimplificationMode), prime_decomposition = false)
        check links.len == 1
        link = links[0]
        check link.crossings.len <= num_crossings
        check link.link_components.len == 3

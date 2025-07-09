import std/unittest

import ../src/links
import ../src/random_links

test "random_link_internal":
    for num_crossings in 1 .. 10:
        var link = random_link_internal[int](num_crossings)
        check link.crossings.len == num_crossings

import std/unittest

import ../src/torus

test "torus_link":
    let T13 = torus_link(1, 3)
    check T13.crossings == []
    let T22 = torus_link(2, 2)
    check T22.crossings == [[5, 4, 7, 6], [1, 0, 3, 2]]
    let T23 = torus_link(2, 3)
    check T23.crossings == [[9, 4, 7, 10], [1, 8, 11, 2], [5, 0, 3, 6]]
    let T33 = torus_link(3, 2)
    check T33.crossings == [[14, 7, 11, 10], [13, 12, 8, 1], [6, 15, 3, 2], [5, 4, 0, 9]]
    # performance:
    # construction of T(10,10): ~10000 per second vs ~300 for Spherogram

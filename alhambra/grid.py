
from typing import Literal, Tuple


def flatgrid_hofromxy(
    x: int, y: int, start_helix: int, start_o: int, p: Literal[9, 10] = 9
) -> Tuple[int, int]:
    if p == 9:
        p = 0
    elif p == 10:
        p = 1
    else:
        raise ValueError
    sx = (p + y) % 2
    sy = (p) % 2
    return (
        start_helix - y + x,
        start_o
        + 23 * (x // 2)
        + 19 * (y // 2)
        + (11 + sx) * (x % 2)
        + (9 + sy) * (y % 2),
    )

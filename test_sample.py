
import pytest
import numpy as np
import sample

def test_add():
    x = 2
    y = 5
    test_sum = add_nums(x, y)
    assert test_sum == 7
import numpy as np
from typing import Tuple


def assert_shape(arr_to_check: np.ndarray, shape: Tuple[int, ...]) -> AssertionError:
    """
    Asserts that the shape of x is equal to shape.
    Parameters
    ----------
    arr_to_check : np.ndarray
        Array to check the shape of.
    shape : Tuple[int, ...]
        Shape to check against.

    Returns
    -------
    AssertionError
        If the shape of x is not equal to shape.

    Examples
    --------
    assert_shape(np.array([1, 2, 3]), (3,))
    """
    assert len(arr_to_check.shape) == len(shape), (arr_to_check.shape, shape)

    for i, (arr_to_check_element, shape_element) in enumerate(zip(arr_to_check.shape, shape)):
        if isinstance(shape_element, int):
            assert arr_to_check_element == shape_element, (f"Index {i} of shape of array to check "
                                                           f"has value {shape[i]}. It should be"
                                                           f" {arr_to_check.shape[i]}")

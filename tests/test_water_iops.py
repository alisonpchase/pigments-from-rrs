import numpy as np
import pytest

from pigments_from_rrs.water_iops import RInw


@pytest.mark.parametrize(
    "kwargs_with_result",
        [
            (dict(lambda_=400, Tc=10, S=30), (1.35, 0.000197)),
            (dict(lambda_=650, Tc=18, S=37), (1.33, 0.000182)),
            (
                dict(lambda_=np.asarray(list(range(400, 701, 50))), Tc=18, S=37),
                (
                    np.asarray([1.35, 1.34, 1.34, 1.34, 1.34, 1.33, 1.33]),
                    np.asarray([0.000193, 0.000189, 0.000187, 0.000185, 0.000183, 0.000182, 0.000180]),
                ),

            ),
        ]
)
def test_RInw(kwargs_with_result):
    def truncate(actual_result):
        truncated_actual_result = tuple(
            [
                (float(str(n)[:4]) if n > 1 else float(str(n)[:8]))
                for n in actual_result
            ]
        )
        return truncated_actual_result

    kwargs, expected_result = kwargs_with_result
    actual_result = RInw(**kwargs)
    if isinstance(kwargs["lambda_"], int):
        assert truncate(actual_result) == expected_result
    else:
        for arr0, arr1 in zip(expected_result, actual_result):
            assert (arr0 == truncate(arr1)).all()

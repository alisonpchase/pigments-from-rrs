import numpy as np
import pytest

from pigments_from_rrs.water_iops import (
    RInw,
    BetaT,
    rho_sw,
    dlnasw_ds,
    PMH,
    betasw124_ZHH2009,
)


def round_result(actual_result):
    truncated_actual_result = tuple(
        [
            (round(n, ndigits=2) if n > 1 else round(n, ndigits=5))
            for n in actual_result
        ]
    )
    return truncated_actual_result


@pytest.mark.parametrize(
    "kwargs_with_result",
        [
            (dict(lambda_=400, Tc=10, S=30), (1.35, 0.00020)),
            (dict(lambda_=650, Tc=18, S=37), (1.34, 0.00018)),
            (
                dict(lambda_=np.asarray(list(range(400, 701, 50))), Tc=18, S=37),
                (
                    np.asarray([1.35, 1.35, 1.34, 1.34, 1.34, 1.34, 1.34]),
                    np.asarray([0.00019, 0.00019, 0.00019, 0.00019, 0.00018, 0.00018, 0.00018]),
                ),

            ),
        ]
)
def test_RInw(kwargs_with_result):
    kwargs, expected_result = kwargs_with_result
    actual_result = RInw(**kwargs)
    if isinstance(kwargs["lambda_"], int):
        assert round_result(actual_result) == expected_result
    else:
        for arr0, arr1 in zip(expected_result, actual_result):
            assert (arr0 == round_result(arr1)).all()


def test_BetaT():
    expected_result = 4.313087388409383e-10
    actual_result = BetaT(Tc=18.0, S=32.0)
    assert expected_result == actual_result


def test_rho_sw():
    expected_result = 1022.9769302697676
    actual_result = rho_sw(Tc=18.0, S=32.0)
    assert expected_result == actual_result


def test_dlnasw_ds():
    expected_result = -0.0005627893517437484
    actual_result = dlnasw_ds(Tc=18.0, S=32.0)
    assert expected_result == actual_result


def test_PMH():
    expected_result = np.asarray(
        [0.90895, 0.87445, 0.87445, 0.87445, 0.87445, 0.84064, 0.84064]
    )
    actual_result = PMH(
        n_wat=np.asarray([1.35, 1.34, 1.34, 1.34, 1.34, 1.33, 1.33])
    )
    assert (expected_result == round_result(actual_result)).all()


def test_betasw124_ZHH2009():
    expected_results = (
        np.asarray([5.2e-04, 3.2e-04, 2.0e-04, 1.4e-04, 0.00009, 0.00007, 0.00005]),
        np.asarray([0.00668, 0.00403, 0.00258, 0.00173, 0.00120, 0.00086, 0.00063]),
        np.asarray([4.1e-04, 2.5e-04, 1.6e-04, 1.1e-04, 7e-05, 5e-05, 4e-05]),
        np.linspace(0.0, 180.0, 18_001),
    )
    actual_results = betasw124_ZHH2009(
        lambda_=np.asarray(list(range(400, 701, 50))),
        S=37,
        Tc=18,
    )
    for arr0, arr1 in zip(expected_results, actual_results):
        if len(arr1) != 18_001:
            assert (arr0 == round_result(arr1)).all()
        else:  # this is theta!
            assert (arr0 == arr1).all()

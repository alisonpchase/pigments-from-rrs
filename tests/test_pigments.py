from pigments_from_rrs.spectra import G1, G2, LNOT, load_data

def test_constants():
    assert G1 == 0.0949
    assert G2 == 0.0794
    assert LNOT == 400


def test_load_data():
    assert "filepath" == load_data(path="filepath")
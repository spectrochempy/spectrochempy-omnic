# ======================================================================================
# Copyright (©) 2015-2025 LCS - Laboratoire Catalyse et Spectrochimie, Caen, France.
# CeCILL-B FREE SOFTWARE LICENSE AGREEMENT
# See full LICENSE agreement in the root directory.
# ======================================================================================
# ruff: noqa

import pytest
from pathlib import Path

from spectrochempy_omnic import OMNICReader, OMNICReaderError


IRDATA = Path(__file__).parent / "data"


def test_read_omnic():
    nd1 = OMNICReader(IRDATA / "nh4y-activation.spg")
    assert nd1.data.shape == (55, 5549)
    assert nd1.x_title == "wavenumbers"
    assert nd1.x_units == "cm^-1"
    assert nd1.y_title == "acquisition timestamp (GMT)"
    assert nd1.y_units == "s"
    assert nd1.units == "absorbance"
    assert nd1.title == "absorbance"
    assert nd1.filename.name == "nh4y-activation.spg"
    assert nd1.y.shape == (55,)
    assert nd1.x.shape == (5549,)
    assert str(nd1) == f"OMNICReader: {nd1.filename.name} {nd1.data.shape}"

    # test read_omnic with byte spg content
    filename_wodger = IRDATA / "wodger.spg"
    with open(filename_wodger, "rb") as fil:
        content = fil.read()
    nd1 = OMNICReader(filename_wodger)
    assert nd1.data.shape == (2, 5549)
    with pytest.raises(OMNICReaderError, match="suffix must be provided in kwargs"):
        # this should raise an error as suffix is not provided
        nd2 = OMNICReader(content)
    nd2 = OMNICReader(content, suffix=".spg")
    assert nd2.data.shape == nd1.data.shape
    assert nd2.x.size == nd2.data.shape[1]

    # Test bytes contents for spa files
    filename = IRDATA / "subdir" / "7_CZ0-100_Pd_101.SPA"
    nds = OMNICReader(filename)
    with open(IRDATA / "subdir" / filename, "rb") as fil:
        content = fil.read()
    nd = OMNICReader(IRDATA / "subdir" / "20-50" / "7_CZ0-100_Pd_21.SPA")
    assert nd.data.shape == (1, 5549)
    nd2 = OMNICReader(IRDATA / "subdir" / "20-50" / "7_CZ0-100_Pd_21.SPA")
    assert nd.data.shape == nds.data.shape == (1, 5549)

    # test import sample IFG
    nd = OMNICReader(
        IRDATA / "carroucell_samp" / "2-BaSO4_0.SPA", interferogram="sample"
    )
    assert nd.interferogram == "sample IFG"
    assert nd.units == "V"
    assert nd.data.shape == (1, 16384)

    # test import background IFG
    nd = OMNICReader(
        IRDATA / "carroucell_samp" / "2-BaSO4_0.SPA", interferogram="background"
    )
    assert nd.interferogram == "background IFG"
    assert nd.units == "V"
    assert nd.data.shape == (1, 16384)

    # import IFG from file without IFG
    a = OMNICReader(
        IRDATA / "subdir" / "20-50" / "7_CZ0-100_Pd_21.SPA", interferogram="sample"
    )
    assert a.data is None

    # rapid_sca series
    a = OMNICReader(IRDATA / "omnic_series/rapid_scan.srs")
    assert a.data.shape == (643, 4160)

    # rapid_sca series, import bg
    a = OMNICReader(IRDATA / "omnic_series/rapid_scan.srs", background=True)
    assert str(a) == f"OMNICReader: {a.filename.name} (1, 4160)"

    # GC Demo
    a = OMNICReader(IRDATA / "omnic_series/GC_Demo.srs")
    assert str(a) == f"OMNICReader: {a.filename.name} (788, 1738)"

    # high speed series
    a = OMNICReader(IRDATA / "omnic_series/high_speed.srs")
    assert str(a) == f"OMNICReader: {a.filename.name} (897, 13898)"

    # high speed series, import bg
    a = OMNICReader(IRDATA / "omnic_series/high_speed.srs", background=True)
    assert str(a) == f"OMNICReader: {a.filename.name} (1, 13898)"

# ======================================================================================
# Copyright (©) 2015-2025 LCS - Laboratoire Catalyse et Spectrochimie, Caen, France.
# CeCILL-B FREE SOFTWARE LICENSE AGREEMENT
# See full LICENSE agreement in the root directory.
# ======================================================================================
# ruff: noqa
import platform
from pathlib import Path

import requests

from spectrochempy_omnic import OMNICReader, OMNICReaderError


def test_read_SOC():
    """upload and read surface optics exemple"""

    baseurl = "https://github.com/chet-j-ski/SOC100_example_data/raw/main/"
    fnames = [
        baseurl + f
        for f in [
            "Fused%20Silica0004.DDR",
            "Fused%20Silica0004.HDR",
            "Fused%20Silica0004.SDR",
        ]
    ]

    ds = OMNICReader(fnames[0])
    assert ds.data.shape == (1, 599)
    assert ds.x_title == "wavenumbers"
    assert ds.x_units == "cm^-1"
    assert (
        str(ds)
        == "OMNICReader: Fused Silica0004.DDR, UnPol, 20.0 [0.0°] 300.0K (1, 599)"
    )

    ds = OMNICReader(fnames[1])
    assert ds.data.shape == (1, 599)
    assert ds.x_title == "wavenumbers"
    assert ds.x_units == "cm^-1"
    assert (
        str(ds)
        == "OMNICReader: Fused Silica0004.HDR, UnPol, 20.0 [0.0°] 300.0K (1, 599)"
    )

    ds = OMNICReader(fnames[2])
    assert ds.data.shape == (1, 599)
    assert ds.x_title == "wavenumbers"
    assert ds.x_units == "cm^-1"
    assert (
        str(ds)
        == "OMNICReader: Fused Silica0004.SDR, UnPol, 20.0 [0.0°] 300.0K (1, 599)"
    )

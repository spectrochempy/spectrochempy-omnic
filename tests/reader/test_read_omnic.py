# ======================================================================================
# Copyright (©) 2015-2025 LCS - Laboratoire Catalyse et Spectrochimie, Caen, France.
# CeCILL-B FREE SOFTWARE LICENSE AGREEMENT
# See full LICENSE agreement in the root directory.
# ======================================================================================
# ruff: noqa

import pytest
from pathlib import Path

from spectrochempy_omnic import OMNICReader, OMNICReaderError


IRDATA = Path(__file__).parent.parent / "data"


class TestOMNICReaderBasics:
    """Tests for basic functionality of OMNICReader."""

    def test_read_spg_file(self):
        """Test reading an SPG file and verify basic properties."""
        spg_file_path = IRDATA / "nh4y-activation.spg"
        nd1 = OMNICReader(spg_file_path)

        # Check data dimensions
        assert nd1.data.shape == (55, 5549)

        # Check metadata
        assert nd1.x_title == "wavenumbers"
        assert nd1.x_units == "cm^-1"
        assert nd1.y_title == "acquisition timestamp (GMT)"
        assert nd1.y_units == "s"
        assert nd1.units == "absorbance"
        assert nd1.title == "absorbance"

        # Check filename and axis shapes
        assert nd1.filename.name == "nh4y-activation.spg"
        assert nd1.y.shape == (55,)
        assert nd1.x.shape == (5549,)

        # Check string representation
        assert str(nd1) == f"OMNICReader: {nd1.filename.name} {nd1.data.shape}"

    def test_read_spg_binary_content(self):
        """Test reading SPG content from binary data."""
        filename_wodger = IRDATA / "wodger.spg"
        with open(filename_wodger, "rb") as fil:
            content = fil.read()

        # Test with file path
        nd1 = OMNICReader(filename_wodger)
        assert nd1.data.shape == (2, 5549)

        # Test with binary content and suffix
        nd2 = OMNICReader(content, suffix=".spg")
        assert nd2.data.shape == nd1.data.shape
        assert nd2.x.size == nd2.data.shape[1]

        # Test error when suffix not provided with binary content
        with pytest.raises(
            OMNICReaderError,
            match="When using bytes content, the suffix must be provided.*",
        ):
            OMNICReader(content)


class TestSPAFiles:
    """Tests for reading SPA files."""

    def test_read_spa_file_different_locations(self):
        """Test reading SPA files from different locations."""
        # File in main subdir
        nds = OMNICReader(IRDATA / "subdir" / "7_CZ0-100_Pd_101.SPA")
        assert nds.data.shape == (1, 5549)

        # File in nested subdir
        nd = OMNICReader(IRDATA / "subdir" / "20-50" / "7_CZ0-100_Pd_21.SPA")
        assert nd.data.shape == (1, 5549)

        # Second read of the same file (tests consistency)
        nd2 = OMNICReader(IRDATA / "subdir" / "20-50" / "7_CZ0-100_Pd_21.SPA")
        assert nd.data.shape == nd2.data.shape


class TestInterferograms:
    """Tests for reading interferogram data from OMNIC files."""

    @pytest.mark.parametrize(
        "ifg_type,expected_units", [("sample", "V"), ("background", "V")]
    )
    def test_read_interferograms(self, ifg_type, expected_units):
        """Test reading different types of interferograms."""
        nd = OMNICReader(
            IRDATA / "carroucell_samp" / "2-BaSO4_0.SPA", interferogram=ifg_type
        )
        assert nd.interferogram == f"{ifg_type} IFG"
        assert nd.units == expected_units
        assert nd.data.shape == (1, 16384)

    def test_read_nonexistent_interferogram(self):
        """Test behavior when trying to read non-existent interferogram."""
        a = OMNICReader(
            IRDATA / "subdir" / "20-50" / "7_CZ0-100_Pd_21.SPA", interferogram="sample"
        )
        assert a.data is None


class TestOMNICSeries:
    """Tests for reading OMNIC series files."""

    @pytest.mark.parametrize(
        "filename,expected_shape,background,expected_bg_shape",
        [
            ("rapid_scan.srs", (643, 4160), True, (1, 4160)),
            ("high_speed.srs", (897, 13898), True, (1, 13898)),
            ("GC_Demo.srs", (788, 1738), False, None),
        ],
    )
    def test_omnic_series(
        self, filename, expected_shape, background, expected_bg_shape
    ):
        """Test reading different types of OMNIC series files."""
        # Test regular series data
        a = OMNICReader(IRDATA / "omnic_series" / filename)
        assert a.data.shape == expected_shape

        # Test background if applicable
        if background:
            bg = OMNICReader(IRDATA / "omnic_series" / filename, background=True)
            assert bg.data.shape == expected_bg_shape
            assert str(bg) == f"OMNICReader: {bg.filename.name} {expected_bg_shape}"


class TestErrorHandling:
    """Tests for error handling in OMNICReader."""

    def test_nonexistent_file(self):
        """Test behavior with non-existent file."""
        nonexistent_path = IRDATA / "nonexistent_file.spg"
        with pytest.raises(OMNICReaderError, match="File not found"):
            OMNICReader(nonexistent_path)

    def test_invalid_suffix(self):
        """Test behavior with invalid file suffix."""
        # Create a temporary text file to test with
        tmp_file = Path(__file__).parent / "temp_test_file.txt"
        try:
            with open(tmp_file, "w") as f:
                f.write("Test content")

            with pytest.raises(OMNICReaderError, match="Invalid suffix"):
                OMNICReader(tmp_file)
        finally:
            # Clean up
            if tmp_file.exists():
                tmp_file.unlink()

    def test_no_suffix(self):
        """Test behavior with no file suffix and no suffix provided."""
        with pytest.raises(OMNICReaderError, match="File has no suffix"):
            # Create a path with no suffix
            path_no_suffix = Path(__file__).parent / "test_file_no_suffix"
            OMNICReader(path_no_suffix)


class TestUtilityFunctions:
    """Tests for utility functions in OMNICReader."""

    def test_is_url(self):
        """Test the is_url function."""
        # This requires accessing the private method, which is not ideal but necessary for coverage
        reader = OMNICReader(IRDATA / "nh4y-activation.spg")
        # Test valid URLs
        assert reader._check_source.__globals__["is_url"]("http://example.com") is True
        assert (
            reader._check_source.__globals__["is_url"]("https://example.com/file.spg")
            is True
        )
        # Test invalid URLs
        assert reader._check_source.__globals__["is_url"]("file.spg") is False
        assert reader._check_source.__globals__["is_url"](123) is False

    def test_utcnow(self):
        """Test the utcnow function."""
        reader = OMNICReader(IRDATA / "nh4y-activation.spg")
        # Just verify it returns a non-None value (actual value will change)
        assert reader._check_source.__globals__["utcnow"]() is not None


class TestReaderMethods:
    """Tests for specific reader methods."""

    def test_history_property(self):
        """Test the history property getter and setter."""
        reader = OMNICReader(IRDATA / "nh4y-activation.spg")

        # Test adding a new history entry
        initial_history_length = len(reader.history)
        reader.history = "Test history entry"
        assert len(reader.history) == initial_history_length + 1
        assert "Test history entry" in reader.history[-1]

        # Test replacing history
        reader.history = []
        assert len(reader.history) == 0

        # Test handling None
        reader.history = None
        assert len(reader.history) == 0

    def test_string_representations(self):
        """Test string representations (__str__ and __repr__)."""
        reader = OMNICReader(IRDATA / "nh4y-activation.spg")

        # Test __str__
        str_representation = str(reader)
        assert "OMNICReader:" in str_representation
        assert reader.filename.name in str_representation
        assert str(reader.data.shape) in str_representation

        # Test __repr__
        repr_representation = repr(reader)
        assert "OMNICReader(" in repr_representation
        assert reader.filename.name in repr_representation
        assert str(reader.data.shape) in repr_representation

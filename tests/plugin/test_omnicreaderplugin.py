# ======================================================================================
# Copyright (©) 2015-2025 LCS - Laboratoire Catalyse et Spectrochimie, Caen, France.
# CeCILL-B FREE SOFTWARE LICENSE AGREEMENT
# See full LICENSE agreement in the root directory.
# ======================================================================================
# ruff: noqa: S101

# Skip tests until main spectrochempy plugin manager is completed
import pytest

try:
    from spectrochempy_omnic.plugin.omnicreaderplugin import OMNICReaderPlugin
except ImportError:
    pytestmark = pytest.mark.skip(
        reason="SpectrochemPy plugin manager not yet finished"
    )


class TestOMNICReaderPlugin:
    def test_get_filetype_info(self):
        """Test that file type info is correctly reported."""
        plugin = OMNICReaderPlugin()
        info = plugin.get_filetype_info()

        assert info["identifier"] == "omnic"
        assert "OMNIC" in info["description"]
        assert "spa" in info["extensions"]
        assert "spg" in info["extensions"]
        assert "srs" in info["extensions"]
        assert info["reader_method"] == "read_omnic"

    def test_read_file_single(self, monkeypatch):
        """Test reading a single file."""
        # Setup mock
        mock_instance = type("MockOMNICInstance", (), {"data": "mock_data"})()

        # Create a mock OMNICReader class
        class MockOMNICReader:
            def __init__(self, file, **kwargs):
                self.file = file
                self.kwargs = kwargs
                self.called_with = {"file": file, **kwargs}

            def __new__(cls, *args, **kwargs):
                return mock_instance

        # Apply the monkeypatch
        monkeypatch.setattr(
            "spectrochempy_omnic.plugin.omnicreaderplugin.OMNICReader", MockOMNICReader
        )

        # Call method
        plugin = OMNICReaderPlugin()
        test_file = "test.spa"
        result = plugin.read_file([test_file])

        # Assertions
        assert len(result) == 1
        assert result[0] is mock_instance

    def test_read_file_multiple(self, monkeypatch):
        """Test reading multiple files."""
        # Setup mocks
        mock_instances = []
        call_args = []

        # Create a mock OMNICReader class
        class MockOMNICReader:
            def __init__(self, file, **kwargs):
                self.file = file
                self.kwargs = kwargs
                self.data = "mock_data"  # Add data attribute
                call_args.append((file, kwargs))
                mock_instances.append(self)

        # Apply the monkeypatch
        monkeypatch.setattr(
            "spectrochempy_omnic.plugin.omnicreaderplugin.OMNICReader", MockOMNICReader
        )

        # Call method
        plugin = OMNICReaderPlugin()
        test_files = ["test1.spa", "test2.spg", "test3.srs"]
        result = plugin.read_file(test_files)

        # Assertions
        assert len(result) == 3
        assert len(mock_instances) == 3
        assert result == mock_instances

        # Verify all files were passed to the reader
        for i, file in enumerate(test_files):
            assert call_args[i][0] == file

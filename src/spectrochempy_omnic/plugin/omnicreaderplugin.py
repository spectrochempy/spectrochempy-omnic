from pluggy import HookimplMarker
from spectrochempy.plugins.readers.readerplugin import ReaderPlugin

from spectrochempy_omnic.reader.read_omnic import OMNICReader

hookimpl = HookimplMarker("spectrochempy")


# generic reader for Omnic files
class OMNICReaderPlugin(ReaderPlugin):
    """Reader for Omnic files."""

    # Hooks implementation
    # --------------------
    @hookimpl
    def get_filetype_info(self):
        return {
            "identifier": "omnic",
            "description": "Nicolet OMNIC files and series (*.spa *.spg *.srs)",
            "extensions": ["spa", "spg", "srs"],
            "reader_method": "read_omnic",
        }

    @hookimpl
    def read_file(self, filenames: list, protocol: str = None, **kwargs) -> object:
        nds = []
        for filename in filenames:
            nd = OMNICReader(filename, protocol=protocol, **kwargs)
            if nd.data is not None:
                nds.append(nd)
        return nds if len(nds) > 0 else None

    @hookimpl
    def install_testdata(self, datadir):
        """Install test data in the specified directory."""
        # from spectrochempy import error_  # here to avoid circular imports
        from spectrochempy_omnic.plugin.utils import download_irdata_directory

        # Download the IR data directory
        download_irdata_directory(datadir)

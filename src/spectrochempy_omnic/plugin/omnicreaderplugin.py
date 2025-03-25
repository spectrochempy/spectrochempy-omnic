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
    def read_file(self, files: list, suffix: str, **kwargs) -> object:
        nds = []
        for file in files:
            nd = OMNICReader(file, **kwargs)
            nds.append(nd)
        return nds

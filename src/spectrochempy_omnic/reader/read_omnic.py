# ======================================================================================
# Copyright (©) 2025 LCS - Laboratoire Catalyse et Spectrochimie, Caen, France.
# CeCILL-B FREE SOFTWARE LICENSE AGREEMENT
# See full LICENSE agreement in the root directory.
# ======================================================================================
# ruff: noqa: S101
"""Module providing a class `OMNICReader` to read OMNIC generated data files."""

__all__ = ["OMNICReader"]

import io
import logging
import re
import struct
import sys
import warnings
from datetime import datetime
from datetime import timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import requests

# ======================================================================
# UTC timezone
# ======================================================================
if sys.version_info >= (3, 11):  # noqa: UP036
    from datetime import UTC
else:
    from datetime import timezone

    UTC = timezone.utc  # noqa: UP017


# ======================================================================
# Exceptions
# ======================================================================
class OMNICReaderError(Exception):
    """Custom exception for OMNIC reader errors."""


class OMNICReaderWarning(Warning):
    """Custom warning for OMNIC reader warnings."""


# ======================================================================
# Utility functions
# ======================================================================
def fromfile(fid, dtype, count):
    """
    Read binary data from a file-like object as a numpy array or scalar.

    This function replaces np.fromfile for io.BytesIO objects.

    Parameters
    ----------
    fid : file-like object
        File-like object to read from.
    dtype : str
        Data type to read ('uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'float32', 'char8').
    count : int
        Number of elements to read.

    Returns
    -------
    numpy.ndarray or scalar
        The data read from the file.
    """
    t = {
        "uint8": "B",
        "int8": "b",
        "uint16": "H",
        "int16": "h",
        "uint32": "I",
        "int32": "i",
        "float32": "f",
        "char8": "c",
    }
    typ = t[dtype] * count
    if dtype.endswith("16"):
        count *= 2
    elif dtype.endswith("32"):
        count *= 4

    out = struct.unpack(typ, fid.read(count))
    if len(out) == 1:
        return out[0]
    return np.array(out)


def is_url(strg):
    """
    Check if a string is a valid URL.

    Parameters
    ----------
    strg : str
        String to check.

    Returns
    -------
    bool
        True if the string is a URL, False otherwise.
    """
    return isinstance(strg, str) and re.match(r"http[s]?:[\/]{2}", strg) is not None


def utcnow():
    """
    Return the current time in UTC with a timezone.

    Returns
    -------
    datetime
        Current UTC time with timezone information.
    """
    if sys.version_info[1] < 12:
        return datetime.utcnow().replace(microsecond=0, tzinfo=ZoneInfo("UTC"))
    return datetime.now(UTC).replace(microsecond=0)


# ======================================================================
# Logger setup
# ======================================================================

# Create logger for the module
logger = logging.getLogger("spectrochempy-omnic")

# Set default level
logger.setLevel(logging.INFO)

# Create console handler if not already attached
if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


# Convenience methods


def info_(msg):
    """
    Log an info message.

    Parameters
    ----------
    msg : str
        Message to log.
    """
    logger.info(msg)


# ======================================================================
# Reader class
# ======================================================================
class OMNICReader:
    """
    Open and convert a Thermo Nicolet OMNIC file.

    Handle Omnic file with these formats: :file:`.spg`, :file:`.spa` or
    :file:`.srs` files

    The reader can also files from Surface Optics Corps, which are compatible with OMNIC.
    Valid extensions are:  :file:`.ddr`, :file:`.hdr` or :file:`.sdr`.

    The collected metadata are:
    - acquisition date (UTC)
    - units of spectrum (absorbance, transmittance, reflectance, Log(1/R),
    Kubelka-Munk, Raman intensity, photoacoustics, volts)
    - units of xaxis (wavenumbers in :math:`cm^{-1}`, wavelengths in nm or micrometer,
    Raman shift in :math:`cm^{-1}`)

    Parameters
    ----------
    source : `str`, `~pathlib.Path` object, byte content object
        The data source to read.
    **kwargs : keyword parameters, optional
        See Other Parameters.

    Other Parameters
    ----------------
    directory : `~pathlib.Path` object objects, optional
        From where to read the files. If not provided or the given path is relative,
        the current working directory is used.
    protocol : `str`, optional
        Protocol parameters is required to read bytes contents as the file type cannot be determined automatically.
        Protocol can be given as e.g ".spg" or "spg" ((Backward compatibility with original reader version in spectrochempy).
    suffix : `str`, optional
        File suffix to read. If not provided, the suffix is determined from the file name.
        (this is equivalent to the protocol parameter). Suffix can be given as e.g ".spg" or "spg".
    interferogram : str, optional
        Apply only to .spa, .hdr, .ddr or .sdr files.
        Default value is None. When set to 'sample' returns the sample interferogram
        of the spa file if present or None if absent. When set to 'background' returns
        the background interferogram of the spa file if present or None if absent.
    return_ifg : str, optional
        Alias of interferogram (Backward compatibility with original reader version in spectrochempy)
    background : bool, optional
        Apply only to .srs files.
        Default value is False. When set to 'True' returns the series background
    return_bg : bool, optional
        Alias of background (Backward compatibility with original reader version in spectrochempy).
    reverse_x : bool, optional
        Apply only to .srs files.
        Specifies whether to reverse the x-axis (wavenumber) data. Default is False.
        In most srs files, the absorbance/intensity data are recorded from high to low
        wavenumbers. However, in some cases the data maybe stored in low to high order.
        In such a case, 'reverse_x' should be set to 'True'.
    """

    suffix = [".spg", ".spa", ".srs", ".ddr", ".hdr", ".sdr"]

    # Initialize attributes for SPG and SPA files
    _timezone = UTC
    description = ""
    data = None
    units = None
    title = None
    name = None
    filename = None
    origin = ""
    original_name = None
    x = None
    x_title = None
    x_units = None
    y = None
    y_timestamp = None
    y_title = None
    y_units = None
    y_labels = None
    mask = None
    interferogram = False
    date = None
    collection_length = None
    collection_length_units = None
    optical_velocity = None
    laser_frequency = None
    laser_frequency_units = None
    _history = []

    # ======================================================================================
    # Initialization and reading
    # ======================================================================================
    def __init__(self, source, **kwargs):
        """
        Initialize the OMNICReader with a data source.

        Parameters
        ----------
        source : str, Path, bytes
            The data source to read.
        **kwargs : dict
            Additional keyword arguments for reading.
        """
        # Check the source
        source, suffix = self._check_source(source, **kwargs)

        # Get the appropriate reader function dynamically
        reader_func = getattr(self, f"_read_{suffix[1:]}")

        # Call the reader function
        reader_func(source, **kwargs)

    @property
    def history(self):
        """
        Get the history of actions made on this class.

        Returns
        -------
        list of str
            List of timestamped history entries.
        """
        history = []
        for date, value in self._history:
            date = date.astimezone(self._timezone).isoformat(
                sep=" ",
                timespec="seconds",
            )
            value = str(value).capitalize()
            history.append(f"{date}> {value}")
        return history

    @history.setter
    def history(self, value):
        """
        Add an entry to the history.

        Parameters
        ----------
        value : str or list
            History entry to add or list of entries to replace current history.
        """
        if value is None:
            return
        if isinstance(value, list):
            # history will be replaced
            self._history = []
            if len(value) == 0:
                return
            value = value[0]
        date = utcnow()
        self._history.append((date, value))

    # ======================================================================================
    # Private methods
    # ======================================================================================
    def __str__(self):
        """
        Return a string representation of the object.

        Returns
        -------
        str
            String representation.
        """
        return f"OMNICReader: {self.filename.name} {self.data.shape}"

    def __repr__(self):
        """
        Return a string representation of the object.

        Returns
        -------
        str
            String representation for developers.
        """
        return f"OMNICReader({self.filename.name}, {self.data.shape})"

    def _read_spg(self, source, **kwargs):
        fid, filename = self._openfid(source, **kwargs)

        # Read name:
        # The name starts at position hex 1e = decimal 30. Its max length
        # is 256 bytes. It is the original filename under which the group has been saved: it
        # won't match with the actual filename if a subsequent renaming has been done in the
        # OS.

        spg_name = self._readbtext(fid, 30, 256)

        # Count the number of spectra
        # From hex 120 = decimal 304, individual spectra are described
        # by blocks of lines starting with "key values",
        # for instance hex[02 6a 6b 69 1b 03 82] -> dec[02 106  107 105 27 03 130]
        # Each of these lines provides positions of data and metadata in the file:
        #
        #     key: hex 02, dec  02: position of spectral header (=> nx, firstx,
        #     lastx, nscans, nbkgscans)
        #     key: hex 03, dec  03: intensity position
        #     key: hex 04, dec  04: user text position
        #     key: hex 1B, dec  27: position of History text
        #     key: hex 64, dec 100: ?
        #     key: hex 66  dec 102: sample interferogram
        #     key: hex 67  dec 103: background interferogram
        #     key: hex 69, dec 105: ?
        #     key: hex 6a, dec 106: ?
        #     key: hex 6b, dec 107: position of spectrum title, the acquisition
        #     date follows at +256(dec)
        #     key: hex 80, dec 128: ?
        #     key: hex 82, dec 130: rotation angle ?
        #
        # the number of line per block may change from file to file but the total
        # number of lines is given at hex 294, hence allowing counting the
        # number of spectra:

        # read total number of lines
        fid.seek(294)
        nlines = fromfile(fid, "uint16", count=1)

        # read "key values"
        pos = 304
        keys = np.zeros(nlines)
        for i in range(nlines):
            fid.seek(pos)
            keys[i] = fromfile(fid, dtype="uint8", count=1)
            pos += 16

        # the number of occurrences of the key '02' is number of spectra
        nspec = np.count_nonzero(keys == 2)

        if nspec == 0:  # pragma: no cover
            raise OMNICReaderError(
                "OSError : File format not recognized - information markers not found",
            )

        # container to hold values
        nx, firstx, lastx = (
            np.zeros(nspec, "int"),
            np.zeros(nspec, "float"),
            np.zeros(nspec, "float"),
        )
        xunits = []
        xtitles = []
        units = []
        titles = []

        # Extracts positions of '02' keys
        key_is_02 = keys == 2  # ex: [T F F F F T F (...) F T ....]'
        indices02 = np.nonzero(key_is_02)  # ex: [1 9 ...]
        position02 = (
            304 * np.ones(len(indices02[0]), dtype="int") + 16 * indices02[0]
        )  # ex: [304 432 ...]

        for i in range(nspec):
            # read the position of the header
            fid.seek(position02[i] + 2)
            pos_header = fromfile(fid, dtype="uint32", count=1)
            # get infos
            info = self._read_header(fid, pos_header, is_first_spectrum=(i == 0))
            nx[i] = info["nx"]
            firstx[i] = info["firstx"]
            lastx[i] = info["lastx"]
            xunits.append(info["xunits"])
            xtitles.append(info["xtitle"])
            units.append(info["units"])
            titles.append(info["title"])

        # check the consistency of xaxis and data units
        if np.ptp(nx) != 0:  # pragma: no cover
            raise OMNICReaderError(
                "Error : Inconsistent data set -"
                " number of wavenumber per spectrum should be "
                "identical",
            )
        if np.ptp(firstx) != 0:  # pragma: no cover
            raise OMNICReaderError(
                "Error : Inconsistent data set - the x axis should start at same value",
            )
        if np.ptp(lastx) != 0:  # pragma: no cover
            raise OMNICReaderError(
                "Error : Inconsistent data set - the x axis should end at same value",
            )
        if len(set(xunits)) != 1:  # pragma: no cover
            raise OMNICReaderError(
                "Error : Inconsistent data set - data units should be identical",
            )
        if len(set(units)) != 1:  # pragma: no cover
            raise OMNICReaderError(
                "Error : Inconsistent data set - x axis units should be identical",
            )
        data = np.ndarray((nspec, nx[0]), dtype="float32")

        # Now the intensity data

        # Extracts positions of '03' keys
        key_is_03 = keys == 3
        indices03 = np.nonzero(key_is_03)
        position03 = 304 * np.ones(len(indices03[0]), dtype="int") + 16 * indices03[0]

        # Read number of spectral intensities
        for i in range(nspec):
            data[i, :] = self._getintensities(fid, position03[i])

        # Get spectra titles & acquisition dates:
        # container to hold values
        spectitles, acquisitiondates, timestamps = [], [], []

        # Extract positions of '6B' keys (spectra titles & acquisition dates)
        key_is_6B = keys == 107
        indices6B = np.nonzero(key_is_6B)
        position6B = 304 * np.ones(len(indices6B[0]), dtype="int") + 16 * indices6B[0]

        # Read spectra titles and acquisition date
        for i in range(nspec):
            # determines the position of informatioon
            fid.seek(position6B[i] + 2)  # go to line and skip 2 bytes
            spa_name_pos = fromfile(fid, "uint32", 1)

            # read omnic filename
            spa_name = self._readbtext(fid, spa_name_pos, 256)
            spectitles.append(spa_name)

            # and the acquisition date
            fid.seek(spa_name_pos + 256)
            timestamp = fromfile(fid, dtype="uint32", count=1)
            # since 31/12/1899, 00:00
            acqdate = datetime(1899, 12, 31, 0, 0, tzinfo=UTC) + timedelta(
                seconds=int(timestamp),
            )
            acquisitiondates.append(acqdate)
            timestamp = acqdate.timestamp()
            # Transform back to timestamp for storage in the Coord object
            # use datetime.fromtimestamp(d, timezone.utc))
            # to transform back to datetime object

            timestamps.append(timestamp)

            # Not used at present
            # -------------------
            # extract positions of '1B' codes (history text), sometimes absent,
            # e.g. peakresolve)
            #  key_is_1B = (keys == 27)
            #  indices1B =  # np.nonzero(key_is_1B)
            #  position1B = 304 * np.ones(len(indices1B[0]), dtype='int') + 16 * indices6B[0]
            #  if len(position1B) != 0:  # read history texts
            #     for j in range(nspec):  determine the position of information
            #        f.seek(position1B[j] + 2)  #
            #        history_pos = fromfile(f,  'uint32', 1)
            #        history =  _readbtext(f, history_pos[0])
            #        allhistories.append(history)

        fid.close()

        # Set object attributes
        self.data = data
        self.units = units[0]
        self.title = titles[0]
        self.name = Path(filename).stem if filename else Path(spg_name).stem
        self.filename = Path(filename) if filename else Path(spg_name)
        self.origin = "omnic"
        self.original_name = spg_name

        # Now get coordinates
        self.x = np.linspace(firstx[0], lastx[0], nx[0])
        self.x_title = xtitles[0]
        self.x_units = xunits[0]

        self.y = np.array(timestamps)  # - min(timestamps)
        self.y_timestamp = min(timestamps)
        self.y_title = "acquisition timestamp (GMT)"
        self.y_units = "s"
        self.y_labels = (acquisitiondates, spectitles)

        self.date = utcnow()
        self.history = f"Imported from spg file {self.filename.name}."

    def _read_spa(self, source, **kwargs):
        fid, filename = self._openfid(source, **kwargs)
        if "return_ifg" in kwargs:
            warnings.warn(
                "The `return_ifg` parameter is deprecated, use `interferogram` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
        interferogram = kwargs.get("interferogram") or kwargs.get("return_ifg")

        # Read name:
        # The name  starts at position hex 1e = decimal 30. Its max length
        # is 256 bytes. It is the original filename under which the spectrum has
        # been saved: it won't match with the actual filename if a subsequent
        # renaming has been done in the OS.
        spa_name = self._readbtext(fid, 30, 256)

        # The acquisition date (GMT) is at hex 128 = decimal 296.
        # Second since 31/12/1899, 00:00
        fid.seek(296)
        timestamp = fromfile(fid, dtype="uint32", count=1)
        acqdate = datetime(1899, 12, 31, 0, 0, tzinfo=UTC) + timedelta(
            seconds=int(timestamp),
        )
        acquisitiondate = acqdate

        # Transform back to timestamp for storage in the Coord object
        # use datetime.fromtimestamp(d, timezone.utc)) to transform back to datetime object
        timestamp = acqdate.timestamp()

        # From hex 120 = decimal 304, the spectrum is described
        # by a block of lines starting with "key values",
        # for instance hex[02 6a 6b 69 1b 03 82] -> dec[02 106  107 105 27 03 130]
        # Each of these lines provides positions of data and metadata in the file:
        #
        #     key: hex 02, dec  02: position of spectral header (=> nx,
        #                                 firstx, lastx, nscans, nbkgscans)
        #     key: hex 03, dec  03: intensity position
        #     #     key: hex 04, dec  04: user text position (custom info, can be present
        #                           several times. The text length is five bytes later)
        #     key: hex 1B, dec  27: position of History text, The text length
        #                           is five bytes later
        #     key: hex 53, dec  83: probably not a position, present when 'Retrieved from library'
        #     key: hex 64, dec 100: ?
        #     key: hex 66  dec 102: sample interferogram
        #     key: hex 67  dec 103: background interferogram
        #     key: hex 69, dec 105: ?
        #     key: hex 6a, dec 106: ?
        #     key: hex 80, dec 128: ?
        #     key: hex 82, dec 130: position of 'Experiment Information', The text length
        #                           is five bytes later. The block gives Experiment filename (at +10)
        #                           Experiment title (+90), custom text (+254), accessory name (+413)
        #     key: hex 92, dec 146: position of 'custom infos', The text length
        #                           is five bytes later.
        #
        # The line preceding the block start with '01' or '0A'
        # The lines after the block generally start with '00', except in few cases where
        # they start by '01'. In such cases, the '53' key is also present
        # (before the '1B').

        # scan "key values"
        pos = 304
        spa_comments = []  # several custom comments can be present
        while "continue":
            fid.seek(pos)
            key = fromfile(fid, dtype="uint8", count=1)

            # print(key, end=' ; ')

            if key == 2:
                # read the position of the header
                fid.seek(pos + 2)
                pos_header = fromfile(fid, dtype="uint32", count=1)
                info = self._read_header(fid, pos_header)

            elif key == 3 and interferogram is None:
                intensities = self._getintensities(fid, pos)

            elif key == 4:
                fid.seek(pos + 2)
                comments_pos = fromfile(fid, "uint32", 1)
                fid.seek(pos + 6)
                comments_len = fromfile(fid, "uint32", 1)
                fid.seek(comments_pos)
                spa_comments.append(fid.read(comments_len).decode("latin-1", "replace"))

            elif key == 27:
                fid.seek(pos + 2)
                history_pos = fromfile(fid, "uint32", 1)
                fid.seek(pos + 6)
                history_len = fromfile(fid, "uint32", 1)
                spa_history = self._readbtext(fid, history_pos, history_len)

            elif key == 102 and interferogram == "sample":
                s_ifg_intensities = self._getintensities(fid, pos)

            elif key == 103 and interferogram == "background":
                b_ifg_intensities = self._getintensities(fid, pos)

            elif key == 00 or key == 1:
                break

            pos += 16

        fid.close()

        if (interferogram == "sample" and "s_ifg_intensities" not in locals()) or (
            interferogram == "background" and "b_ifg_intensities" not in locals()
        ):
            info_("No interferogram found, read_spa returns None")
            return
        if interferogram == "sample":
            intensities = s_ifg_intensities
        elif interferogram == "background":
            intensities = b_ifg_intensities

        # load intensity into the  NDDataset
        self.data = np.array(intensities[np.newaxis], dtype="float32")

        if interferogram == "background":
            title = "sample acquisition timestamp (GMT)"  # bckg acquisition date is not known for the moment...
        else:
            title = "acquisition timestamp (GMT)"  # no ambiguity here

        self.y = np.array([timestamp], dtype="float32")
        self.y_timestamp = timestamp
        self.y_title = title
        self.y_units = "s"
        self.y_labels = ([acquisitiondate], [spa_name])

        # useful when a part of the spectrum/ifg has been blanked:
        self.mask = np.isnan(self.data)

        self.interferogram = interferogram is not None
        self.name = Path(filename).stem if filename else Path(spa_name).stem
        self.filename = Path(filename) if filename else Path(spa_name)
        self.origin = "omnic"
        self.original_name = spa_name

        if interferogram is None:
            self.interferogram = False
            self.units = info["units"]
            self.title = info["title"]

            # now add coordinates
            nx = info["nx"]
            firstx = info["firstx"]
            lastx = info["lastx"]
            xunit = info["xunits"]
            xtitle = info["xtitle"]

            self.x = np.linspace(firstx, lastx, int(nx))
            self.x_title = xtitle
            self.x_units = xunit

        else:  # interferogram
            if interferogram == "sample":
                self.interferogram = "sample IFG"
            else:
                self.interferogram = "background IFG"
            self.units = "V"
            self.title = "detector signal"

            self.x = np.arange(len(intensities))
            self.x_title = "data points"
            self.x_units = None

        self.name = self.original_name  # to be consistent with omnic behaviour

        self.collection_length = info["collection_length"] / 100
        self.collection_length_units = "s"
        self.optical_velocity = info["optical_velocity"]
        self.laser_frequency = info["reference_frequency"]
        self.laser_frequency_units = "cm^-1"

        if len(spa_comments) > 1:
            self.description += "# Comments from Omnic:\n"
            for comment in spa_comments:
                self.description += comment + "\n---------------------\n"

        self.history = "Imported from spa file(s)"

        if "spa_history" in locals() and len("spa_history".strip(" ")) > 0:
            self.history = (
                "Data processing history from Omnic :\n------------------------------------\n"
                + spa_history
            )

        self.date = utcnow()

    def _read_srs(self, source, **kwargs):
        fid, filename = self._openfid(source, **kwargs)
        if "return_bg" in kwargs:
            warnings.warn(
                "The `return_bg` parameter is deprecated, use `background` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
        background = kwargs.get("background", False) or kwargs.get("return_bg", False)
        reverse_x = kwargs.get("reverse_x", False)

        # read the file and determine whether it is a rapidscan or a high speed real time
        is_rapidscan, is_highspeed, is_tg = False, False, False

        """ At pos=304 (hex:130) is the position of the '02' key for series. Here we don't use it.
        Instead, we use one of the following sequence :

        RapidScan series:
        ----------------
        the following sequence appears 3 times in the file
        b'\x02\x00\x00\x00\x18\x00\x00\x00\x00\x00\x48\x43\x00\x50\x43\x47

        They are used to assert the srs file is rapid_scan and to locate headers and data:
        - The 1st one is located 152 bytes after the series header position
        - The 2nd one is located 152 bytes before the background header position and
        56 bytes before either the background data / or the background title and infos
        followed by the background data
        - The 3rd one is located 60 bytes before the series data (spectre/ifg names and
        intensities


        High Speed Real time series:
        ---------------------------
        the following sequence appears 4 times in the file:
        b"\x02\x00\x00\x00\x18\x00\x00\x00\x00\x00\x48\x43\x00\xc8\xaf\x47"

        They are used to assert the srs file is high speed and to locate headers and data:
        - The 1st one is located 152 bytes after the series header position
        - The 2nd one is located 152 bytes before the background header position and
        56 bytes before either the background data / or the background title and infos
        followed by the background data
        - The 3rd one is located 60 bytes before some data (don't know yet what it is)
        - The 4th one is located 60 bytes before the series data (spectra)

        TGA/IR or GC series:
        ---------------------------
        the following sequence appears 3 times in the file:
        b"\x02\x00\x00\x00\x18\x00\x00\x00\x00\x00", the next bytes can differ from
        one file to another.

        As it is common to other types of TG/IR they will be used to assert if the srs file
        is TGA/IR or GC  *after the other formats*. They also allows locating headers and
        data:
        - The 1st one is located 152 bytes after the series header position
        - The 2nd one is located 152 bytes before the background header position and
        56 bytes before either the background data / or the background title and infos
        followed by the background data
        - The 3rd one is located 60 bytes before the series data (spectre/ifg names and
        intensities ?
        """

        sub_rs = b"\x02\x00\x00\x00\x18\x00\x00\x00\x00\x00\x48\x43\x00\x50\x43\x47"
        sub_hs = b"\x02\x00\x00\x00\x18\x00\x00\x00\x00\x00\x48\x43\x00\xc8\xaf\x47"
        sub_tg = b"\x02\x00\x00\x00\x18\x00\x00\x00\x00\x00"

        # find the first occurence and determine whether the srs is rapidscan or high
        # speed real time
        fid.seek(0)
        bytestring = fid.read()

        # try rapidscan first:
        pos = bytestring.find(sub_rs, 1)
        if pos > 0:
            is_rapidscan = True
        else:
            # not rapidscan, try high speed real time
            pos = bytestring.find(sub_hs, 1)
            if pos > 0:
                is_highspeed = True
            else:
                # neith rapid scan nor high speed real time, try TGA/IR
                pos = bytestring.find(sub_tg, 1)
                if pos > 0:
                    is_tg = True

                else:
                    raise OMNICReaderError(
                        "The reader is only implemented for Rapid Scan, "
                        "High Speed Real Time, GC or TGA srs files. If you think "
                        "your file belongs to one of these types, or if "
                        "you'd like an update of the reader to read your "
                        "file type, please report the issue on "
                        "https://github.com/spectrochempy/spectrochempy-omnic"
                        "/issues ",
                    )

        if is_rapidscan:
            # determine whether the srs is reprocessed. At pos=292 (hex:124) appears a
            # difference between pristine and reprocessed series
            fid.seek(292)
            key = fromfile(fid, dtype="uint8", count=16)[0]
            if key == 39:  # (hex: 27)
                is_reprocessed = False
            elif key == 15:  # (hex = 0F)
                is_reprocessed = True
            else:
                raise OMNICReaderError(
                    "The file is not recognized as a Rapid Scan "
                    "srs file. Please report the issue on "
                    "https://github.com/spectrochempy/spectrochempy-omnic"
                    "/issues ",
                )

            # find the 2 following starting indexes of sub_rs.
            # we will use the 1st (-> series info), the 2nd (-> background) and
            # the 3rd  (-> data)

            fid.seek(0)
            bytestring = fid.read()
            index = [pos]
            while pos != -1:
                pos = bytestring.find(sub_rs, pos + 1)
                index.append(pos)

            index = np.array(index[:-1]) + [-152, -152, 60]

            if len(index) != 3:
                raise OMNICReaderError(
                    "The file is not recognized as a Rapid Scan "
                    "srs file. Please report the issue on "
                    "https://github.com/spectrochempy/spectrochempy-omnic"
                    "/issues ",
                )

            pos_info_data = index[0]
            pos_info_bg = index[1]
            pos_data = index[2]

            # read series data, except if the user asks for the background
            if not background:
                info = self._read_header(fid, pos_info_data)
                names, data = self._read_srs_spectra(
                    fid, pos_data, info["ny"], info["nx"]
                )

                # now get series history
                if not is_reprocessed:
                    history = info["history"]
                else:
                    # In reprocessed series the updated "DATA PROCESSING HISTORY" is located right after
                    # the following 16 byte sequence:
                    sub = b"\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff"
                    pos = bytestring.find(sub) + 16
                    history = self._readbtext(fid, pos, None)

            # read the background if the user asked for it.
            if background:
                # First get background info
                info = self._read_header(fid, pos_info_bg)

                if "background_name" not in info:
                    # it is a short header
                    fid.seek(index[1] + 208)
                    data = fromfile(fid, dtype="float32", count=info["nx"])
                else:
                    # longer header, in such case the header indicates a spectrum
                    # but the data are those of an ifg... For now need more examples
                    return

                # uncomment below to load the last datafield has the same dimension as the time axis
                # its function is not known. related to Grams-schmidt ?

                # pos = _nextline(pos)
                # found = False
                # while not found:
                #     pos += 16
                #     f.seek(pos)
                #     key = fromfile(f, dtype='uint8', count=1)
                #     if key == 1:
                #         pos += 4
                #         f.seek(pos)
                #         X = fromfile(f, dtype='float32', count=info['ny'])
                #         found = True

        if is_highspeed:
            # find the 3 following starting indexes of sub.
            # 1st -> series info),
            # 2nd -> background ?
            # 3rd -> data ?
            # 4th  -> ?
            fid.seek(0)
            bytestring = fid.read()
            index = [pos]
            while pos != -1:
                pos = bytestring.find(sub_hs, pos + 1)
                index.append(pos)

            index = np.array(index[:-1]) + [-152, -152, 0, 60]

            pos_info_data = index[0]
            pos_bg = index[1]
            pos_x = index[2]
            pos_data = index[3]

            if len(index) != 4:
                raise OMNICReaderError(
                    "The file is not recognized as a High Speed Real "
                    "Time srs file. Please report the issue on "
                    "https://github.com/spectrochempy/spectrochempy-omnic"
                    "/issues ",
                )

            if not background:
                info = self._read_header(fid, pos_info_data)
                # container for names and data

                names, data = self._read_srs_spectra(
                    fid, pos_data, info["ny"], info["nx"]
                )

                # Get series history. on the sample file, the history seems overwritten by
                # some post-processing, so info["history"] returns a corrupted string.
                # The "DATA PROCESSING HISTORY" (as indicated by omnic) is located right
                # after the following 16 byte sequence:
                sub = (
                    b"\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff"
                )
                pos = bytestring.find(sub) + 16
                history = self._readbtext(fid, pos, None)

                # read the background if the user asked for it.

            elif background:
                # First get background info
                info = self._read_header(fid, pos_bg)

                if "background_name" not in info:
                    # it is a short header
                    fid.seek(index[1] + 208)
                    data = fromfile(fid, dtype="float32", count=info["nx"])
                else:
                    # longer header, in such case the header indicates a spectrum
                    # but the data are those of an ifg... For now need more examples
                    return

        if is_tg:
            fid.seek(0)
            bytestring = fid.read()
            index = [pos]
            while pos != -1:
                pos = bytestring.find(sub_tg, pos + 1)
                index.append(pos)

            index = np.array(index[:-1]) + [-152, -152, 60]

            if len(index) != 3:
                raise OMNICReaderError(
                    "The file is not recognized as a TG IR or GC "
                    "srs file. Please report the issue on "
                    "https://github.com/spectrochempy/spectrochempy-omnic"
                    "/issues ",
                )

            pos_info_data = index[0]
            pos_info_bg = index[1]
            pos_data = index[2]

            # read series data, except if the user asks for the background
            if not background:
                info = self._read_header(fid, pos_info_data)
                names, data = self._read_srs_spectra(
                    fid, pos_data, info["ny"], info["nx"]
                )
                # Note: info["history"] is empty in TG IR or GC series
                # the position of the history is indiated at pos 856 or 878 depending on the
                # file.

            # read the background if the user asked for it.
            if background:
                # First get background info
                info = self._read_header(fid, pos_info_bg)

                if "background_name" not in info:
                    # it is a short header
                    fid.seek(index[1] + 208)
                    data = fromfile(fid, dtype="float32", count=info["nx"])
                else:
                    # longer header, in such case the header indicates a spectrum
                    # but the data are those of an ifg... For now need more examples
                    return

        if not background:
            self.data = data if not reverse_x else data[:, ::-1]
        else:
            self.data = np.expand_dims(data, axis=0)

        # in case part of the spectra/ifg has been blanked:
        self.mask = np.isnan(self.data)

        self.units = info["units"]
        self.title = info["title"]
        self.name = Path(filename).stem if filename else "unnamed"
        self.filename = Path(filename) if filename else "unknown"
        self.origin = "omnic"

        # now add coordinates
        self.x = np.linspace(info["firstx"], info["lastx"], int(info["nx"]))
        self.x_title = info["xtitle"]
        self.x_units = info["xunits"]

        # specific infos for series data
        if not background:
            self.name = info["name"]
            self.y = np.around(
                np.linspace(info["firsty"], info["lasty"], info["ny"]), 3
            )
            self.y_title = "Time"
            self.y_units = "minute"
            self.y_labels = names

        else:
            self.y = np.array([0], dtype="float32")
        self.y_timestamp = min(self.y)
        if "history" in locals():
            self.history = (
                "Omnic 'DATA PROCESSING HISTORY' :\n"
                "--------------------------------\n" + history,
            )
        self.history = " imported from srs file " + str(self.filename)

        self.laser_frequency = info["reference_frequency"]
        self.laser_frequency_units = "cm^-1"
        self.collection_length = info["collection_length"]
        self.collection_length_units = "s"

        if self.x_units is None and self.x_title == "data points":
            # interferogram
            self.interferogram = True

        fid.close()

    def _read_ddr(self, *args, **kwargs):
        self._read_spa(*args, **kwargs)
        self.history[-1] = "Imported from ddr file(s)"

    def _read_hdr(self, *args, **kwargs):
        self._read_spa(*args, **kwargs)
        self.history[-1] = "Imported from hdr file(s)"

    def _read_sdr(self, *args, **kwargs):
        self._read_spa(*args, **kwargs)
        self.history[-1] = "Imported from sdr file(s)"

    def _get_suffix_from_kwargs(self, **kwargs):
        suffix = kwargs.get("protocol") or kwargs.get("suffix")
        if suffix is None:
            return None
        if not suffix.startswith("."):
            suffix = "." + suffix
        return suffix.lower()

    def _check_source(self, source, **kwargs):
        """
        Check and validate the source type and suffix.

        Parameters
        ----------
        source : str, Path, bytes
            The data source to check.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        tuple
            Processed source and validated suffix.

        Raises
        ------
        OMNICReaderError
            If the source or suffix is invalid.
        """
        # Check the source type and validate the suffix if necessary
        kw_suffix = self._get_suffix_from_kwargs(**kwargs)

        if is_url(source):
            suffix = Path(source).suffix.lower()
            return source, kw_suffix if suffix not in self.suffix else suffix

        if not isinstance(source, bytes | dict):
            # Check if source is a string or Path object
            if not isinstance(source, Path):
                try:
                    source = Path(source)
                except Exception as e:
                    raise OMNICReaderError(
                        f"Invalid type for path: {type(source)}. Expected str or Pathlib object"
                    ) from e

            # Check if the file has a valid suffix
            if not source.suffix and kw_suffix is None:
                raise OMNICReaderError(
                    f"File has no suffix and suffix (or protocol) parameter is not provided in kwargs. Expected one of {self.suffix}"
                )

            # Check if source is a valid file path
            if not source.exists():
                raise OMNICReaderError(f"File not found: {source}")

            suffix = source.suffix.lower() if source.suffix else kw_suffix
            if suffix not in self.suffix:
                raise OMNICReaderError(
                    f"Invalid suffix: {suffix}. Expected one of {self.suffix}"
                )

        else:  # source is a content
            # Check if suffix is provided
            suffix = kw_suffix
            if isinstance(source, dict) and suffix is None:
                # dict must content a single element
                if len(source) != 1:
                    raise OMNICReaderError(
                        "When using dict content, the dict must contain a single element. (multiple content reading not yet supported)"
                    )
                suffix = Path(list(source.keys())[0]).suffix.lower()
            if suffix is None:
                raise OMNICReaderError(
                    "When using bytes content, the suffix must be provided as a filename suffix or "
                    "in kwargs (using protocol or suffix parameter)."
                )
            # Validate the provided suffix
            if suffix not in self.suffix:
                raise OMNICReaderError(
                    f"Invalid suffix: {suffix}. Expected one of {self.suffix}"
                )

        return source, suffix

    def _openfid(self, source, mode="rb", **kwargs):
        """
        Open a file-like object from various source types.

        Parameters
        ----------
        source : str, Path, bytes, dict
            The data source to open.
        mode : str, optional
            The file opening mode, default is 'rb'.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        file-like object
            Opened file object.
        """
        # default encoding
        encoding = "utf-8"

        content = None
        filename = None

        if is_url(source):
            # use request to read the remote content
            r = requests.get(source, allow_redirects=True, timeout=10)
            r.raise_for_status()
            content = r.content
            encoding = r.encoding
            filename = Path(source).name

        elif isinstance(source, bytes):
            # use BytesIO to read the bytes content
            content = source

        elif isinstance(source, dict):
            content = list(source.values())[0]
            filename = list(source.keys())[0]

        else:
            # transform filename to a Path object if not yet the case
            filename = Path(source)

        # Create the file ID
        if content is not None:
            # if a content has been passed
            fid = (
                io.BytesIO(content)
                if mode == "rb"
                else io.StringIO(content.decode(encoding))
            )
        else:
            fid = open(filename, mode=mode)  # noqa: SIM115

        return fid, filename

    def _read_header(self, fid, pos,is_first_spectrum=True):
        r"""
        Read spectrum/ifg/series header.

        Parameters
        ----------
        fid : BufferedReader
            The buffered binary stream.

        pos : int
            The position of the header (see Notes).

        is_first_spectrum : bool, optional
            Indicates if this is the first spectrum being read. Default is True.

        Returns
        -------
            dict, int
            Dictionary and current position in file

        Notes
        -----
            So far, the header structure is as follows:

            - starts with b'\x01' , b'\x02', b'\x03' ... maybe indicating the header "type"
            - nx (UInt32): 4 bytes behind
            - xunits (UInt8): 8 bytes behind. So far, we have the following correspondence:

                * `x\01` : wavenumbers, cm-1
                * `x\02` : datapoints (interferogram)
                * `x\03` : wavelength, nm
                * `x\04' : wavelength, um
                * `x\20' : Raman shift, cm-1

            - data units (UInt8): 12 bytes behind. So far, we have the following
            correspondence:

                * `x\0B` : reflectance (%)
                * `x\0C` : Kubelka_Munk
                * `x\0F` : single beam
                * `x\11` : absorbance
                * `x\10` : transmittance (%)
                * `x\16` : Volts (interferogram)
                * `x\1A` : photoacoustic
                * `x\1F` : Raman intensity

            - first x value (float32), 16 bytes behind
            - last x value (float32), 20 bytes behind
            - ... unknown
            - scan points (UInt32), 28 bytes behind
            - zpd (UInt32),  32 bytes behind
            - number of scans (UInt32), 36 bytes behind
            - ... unknown
            - number of background scans (UInt32), 52 bytes behind
            - ... unknown
            - collection length in 1/100th of sec (UIint32), 68 bytes behind
            - ... unknown
            - reference frequency (float32), 80 bytes behind
            - ...
            - optical velocity (float32), 188 bytes behind
            - ...
            - spectrum history (text), 208 bytes behind

            For "rapid-scan" srs files:

            - series name (text), 938 bytes behind
            - collection length (float32), 1002 bytes behind
            - last y (float 32), 1006 bytes behind
            - first y (float 32), 1010 bytes behind
            - ny (UInt32), 1026
            - ... y unit could be at pos+1030 with 01 = minutes ?
            - history (text), 1200 bytes behind (only initial history.
            When reprocessed, updated history is at the end of the file after the
            b`\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff` sequence

        """
        out = {}
        # determine the type of file
        fid.seek(0)
        bytes = fid.read(18)
        if bytes == b"Spectral Data File":
            filetype = "spa, spg"
        elif bytes == b"Spectral Exte File":
            filetype = "srs"

        # nx
        fid.seek(pos + 4)
        out["nx"] = fromfile(fid, "uint32", count=1)

        # xunits
        fid.seek(pos + 8)
        key = fromfile(fid, dtype="uint8", count=1)
        if key == 1:
            out["xunits"] = "cm^-1"
            out["xtitle"] = "wavenumbers"
        elif key == 2:
            out["xunits"] = None
            out["xtitle"] = "data points"
        elif key == 3:  # pragma: no cover
            out["xunits"] = "nm"
            out["xtitle"] = "wavelengths"
        elif key == 4:  # pragma: no cover
            out["xunits"] = "um"
            out["xtitle"] = "wavelengths"
        elif key == 32:  # pragma: no cover
            out["xunits"] = "cm^-1"
            out["xtitle"] = "raman shift"
        else:  # pragma: no cover
            out["xunits"] = None
            out["xtitle"] = "xaxis"
            info_("The nature of x data is not recognized, xtitle is set to 'xaxis'")

        # data units
        fid.seek(pos + 12)
        key = fromfile(fid, dtype="uint8", count=1)
        if key == 17:
            out["units"] = "absorbance"
            out["title"] = "absorbance"
        elif key == 16:  # pragma: no cover
            out["units"] = "percent"
            out["title"] = "transmittance"
        elif key == 11:  # pragma: no cover
            out["units"] = "percent"
            out["title"] = "reflectance"
        elif key == 12:  # pragma: no cover
            out["units"] = None
            out["title"] = "log(1/R)"
        elif key == 15:  # pragma: no cover
            out["units"] = None
            out["title"] = "single beam"
        elif key == 20:  # pragma: no cover
            out["units"] = "Kubelka_Munk"
            out["title"] = "Kubelka-Munk"
        elif key == 21:
            out["units"] = None
            out["title"] = "reflectance"
        elif key == 22:
            out["units"] = "V"
            out["title"] = "detector signal"
        elif key == 26:  # pragma: no cover
            out["units"] = None
            out["title"] = "photoacoustic"
        elif key == 31:  # pragma: no cover
            out["units"] = None
            out["title"] = "Raman intensity"
        else:  # pragma: no cover
            out["units"] = None
            out["title"] = "intensity"
            if is_first_spectrum:
                info_(f"The nature of data is not recognized (key == {key}), title set to 'Intensity'")

        # firstx, lastx
        fid.seek(pos + 16)
        out["firstx"] = fromfile(fid, "float32", 1)
        fid.seek(pos + 20)
        out["lastx"] = fromfile(fid, "float32", 1)
        fid.seek(pos + 28)

        out["scan_pts"] = fromfile(fid, "uint32", 1)
        fid.seek(pos + 32)
        out["zpd"] = fromfile(fid, "uint32", 1)
        fid.seek(pos + 36)
        out["nscan"] = fromfile(fid, "uint32", 1)
        fid.seek(pos + 52)
        out["nbkgscan"] = fromfile(fid, "uint32", 1)
        fid.seek(pos + 68)
        out["collection_length"] = fromfile(fid, "uint32", 1)
        fid.seek(pos + 80)
        out["reference_frequency"] = fromfile(fid, "float32", 1)
        fid.seek(pos + 188)
        out["optical_velocity"] = fromfile(fid, "float32", 1)

        if filetype == "spa, spg":
            out["history"] = self._readbtext(fid, pos + 208, None)

        if filetype == "srs":
            if out["nbkgscan"] == 0 and out["firstx"] > out["lastx"]:
                # an interferogram in rapid scan mode
                out["firstx"], out["lastx"] = out["lastx"], out["firstx"]

            out["name"] = self._readbtext(fid, pos + 938, 256)
            # Hack because name seems not to be well read for srs
            out["name"] = out["name"].split("\n")[0]
            fid.seek(pos + 1002)
            out["collection_length"] = fromfile(fid, "float32", 1) * 60
            fid.seek(pos + 1006)
            out["lasty"] = fromfile(fid, "float32", 1)
            fid.seek(pos + 1010)
            out["firsty"] = fromfile(fid, "float32", 1)
            fid.seek(pos + 1026)
            out["ny"] = fromfile(fid, "uint32", 1)
            #  y unit could be at pos+1030 with 01 = minutes ?
            out["history"] = self._readbtext(fid, pos + 1200, None)

            if self._readbtext(fid, pos + 208, 256)[:10] == "Background":
                # it is the header of a background
                out["background_name"] = self._readbtext(fid, pos + 208, 256)[10:]

        return out

    def _read_srs_spectra(self, fid, pos_data, n_spectra, n_points):
        """
        Read the spectra/interferogram names and data of a series.

        fid: BufferedReader
        pos_data: int
        n_spectra: int
        n_points: int

        returns: names (list), spectral data (ndarray)
        """
        # container for names and data
        names = []
        data = np.zeros((n_spectra, n_points))

        # read the spectra/interferogram names and data
        # the first one....
        pos = pos_data
        names.append(self._readbtext(fid, pos, 256))
        pos += 84
        fid.seek(pos)
        data[0, :] = fromfile(fid, dtype="float32", count=n_points)[:]
        pos += n_points * 4
        # ... and the remaining ones:
        for i in np.arange(n_spectra)[1:]:
            pos += 16
            names.append(self._readbtext(fid, pos, 256))
            pos += 84
            fid.seek(pos)
            data[i, :] = fromfile(fid, dtype="float32", count=n_points)[:]
            pos += n_points * 4

        return names, data

    @staticmethod
    def _readbtext(fid, pos, size):
        # Read some text in binary file of given size. If size is None, the etxt is read
        # until b\0\ is encountered.
        # Returns utf-8 string
        fid.seek(pos)
        if size is None:
            btext = b""
            while fid.read(1) != b"\x00":
                btext += fid.read(1)
        else:
            btext = fid.read(size)
        btext = re.sub(b"\x00+", b"\n", btext)

        if btext[:1] == b"\n":
            btext = btext[1:]

        if btext[-1:] == b"\n":
            btext = btext[:-1]

        try:
            text = btext.decode(encoding="utf-8")  # decode btext to string
        except UnicodeDecodeError:
            try:
                text = btext.decode(encoding="latin_1")
            except UnicodeDecodeError:  # pragma: no cover
                text = btext.decode(encoding="utf-8", errors="ignore")
        return text

    @staticmethod
    def _nextline(pos):
        # reset current position to the beginning of next line (16 bytes length)
        return 16 * (1 + pos // 16)

    @staticmethod
    def _getintensities(fid, pos):
        # get intensities from the 03 (spectrum)
        # or 66 (sample ifg) or 67 (bg ifg) key,
        # returns a ndarray

        fid.seek(pos + 2)  # skip 2 bytes
        intensity_pos = fromfile(fid, "uint32", 1)
        fid.seek(pos + 6)
        intensity_size = fromfile(fid, "uint32", 1)
        nintensities = int(intensity_size / 4)

        # Read and return spectral intensities
        fid.seek(intensity_pos)
        return fromfile(fid, "float32", int(nintensities))


if __name__ == "__main__":
    pass

# SpectroChemPy OMNIC Reader

A standalone reader (and soon a plugin for SpectroChemPy) that enables
reading Thermo Scientific™ OMNIC™ spectroscopy files (*.spa,*.spg, *.srs) and Surface Optics Corp. Files (*.hdr, *.sdr,*.ddr)

## Installation

```bash
pip install spectrochempy-omnic
```

## Usage

```python
from spectrochempy_omnic import OMNICReader as read

# Read an OMNIC file
res = read("path/to/your/file.spg")

# get results numpy arrays
data = res.data  # array of data
shape = data.shape

# get axis
x = res.x
y = res.y

xunits, xtitle = res.x_units, res.x_title
yunits, ytitle = res.y_units, res.y_title

```

## Requirements

- Python >=3.10
- NumPy

## License

This project is licensed under the CeCILL-B FREE SOFTWARE LICENSE AGREEMENT.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

If you encounter any problems, please [file an issue](https://github.com/spectrochempy/spectrochempy-omnic/issues) along with a detailed description.

## Authors

- Arnaud Travert (<contact@spectrochempy.fr>)
- Christian Fernandez (<contact@spectrochempy.fr>)

## Citation

If you use this software in your research, please cite:

```bibtex
@software{spectrochempy_omnic,
    title = {SpectroChemPy OMNIC Reader},
    author = {Travert, Arnaud and Fernandez, Christian},
    url = {https://github.com/spectrochempy/spectrochempy-omnic},
    version = {0.1.0},
    year = {2025}
}

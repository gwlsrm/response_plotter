# Detector response polynomial plotter
## About
Script plots detector response polynomials from LSRM EffMaker detectors database (http://lsrm.ru/en/products/detail.php?ELEMENT_CODE=EffMaker)

## How to use
Run `DRGen_P.exe` from EffMaker package. Go to menu `File`->`Database` (ctrl+B) and export response from selected detector to csv-file.
Run `plot_detector_response_poly` from command-line with parameters:

- this shows available energies, for which response was calculated:
```sh
python plot_detector_response_poly.py <filename.csv>
```
- this plots response polynomial for selected energy
```bash
python plot_detector_response_poly.py <filename.csv> <point_num>
```

## Dependencies 
- `matplotlib`
- `numpy`


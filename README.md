# Closed-loop Resonant Axial-Scanning High-speed 2-photon (CRASH2p) microscopy data analysis

This repository contains data analysis codes used in the paper [McNulty, Wu, et. al., 2025](https://doi.org/10.1038/s41467-025-60648-x).

# Getting Started

### Software requirements

[MATLAB](https://www.mathworks.com/products/matlab.html) R2024a or later required, with the following toolboxes: [Signal Processing](https://www.mathworks.com/products/signal.html), [Curve Fitting](https://www.mathworks.com/products/curvefitting.html), [Statistics and Machine Learning](https://www.mathworks.com/products/statistics.html), [Optimization](https://www.mathworks.com/products/optimization.html), and [Image Processing](https://www.mathworks.com/products/image-processing.html).

[SLEAP](https://sleap.ai/) required for behavioral image labeling, which can be done separately and does not require this code repository, as long as the labels are saved with the correct file path pattern; see `demo_hyperscope_pipeline` for details.

### Sample data and demo

Start with the script `demo_hyperscope_pipeline` under the root folder. Sample datasets can be found in the Dataverse repository accompanying our paper: https://doi.org/10.7910/DVN/ZNJ8U9.

# Citation

If you use codes provided in this repository, please cite:

> McNulty, P., Wu, R., Yamaguchi, A., Heckscher, E. S., Haas, A., Nwankpa, A., ... & Gershow, A. M. (2025). Closed-loop two-photon functional imaging in a freely moving animal. Nature Communications, 16(1), 5950. [[paper]](https://doi.org/10.1038/s41467-025-60648-x)
> 

# License

Codes provided in this repository are free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

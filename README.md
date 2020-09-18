# FLOWPortrait

Code to perform FLOW Portrait based analysis and to support paper by N.Linden, D.Tabuena, N.Steinmetz, S.Brunton, B.Brunton.

## Software Requirements
- Developed using MATLAB R2017b. 
- Requires MATLAB's Image Processing Toolbox (<https://www.mathworks.com/products/image.html>)

## License

This code is available under the MIT License - see LICENSE.md for details

## Acknowledgements

- Kristjan Onu for developing LCS Tool, a software to compute Finite Time Lyapunov Exponents [K. Onu et al.](https://www.sciencedirect.com/science/article/abs/pii/S187775031400163X#!)
- Mohd Kharbat for developing Horn-Schunck Optical Flow Method, a MATLAB implementation of the Horn-Schunck optical flow method [Mohd Kharbat](https://www.mathworks.com/matlabcentral/fileexchange/22756-horn-schunck-optical-flow-method)

## Code Organization

The software is organized into three directories (**flow_portraits/**, **figures/**, and **data/** ) as follows:

- The directory **flow_portraits/** contains all code necesary to perform a FLOW portrait based analysis.

  - **flow_portraits/** contains additional required code 
  - **flow_portraits/LCS-tool/** contains required code from the LCS Tool (see Acknowledgements above for details) to compute the FTLE
  - **flow_portraits/demo/** which outlines steps to compute a FLOW portrait and to choose an integration length

- The directory **figures/** contains scripts to compute the panes in key figures for the FLOW portrait paper

- The directory **data/** contains all data necesary to run the demo scripts and to compute the key figures

  


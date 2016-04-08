# compressed-sensing-diffusion-lung-MRI

This repository contains data, code and results for a , presented in the paper **Incorporation of prior knowledge of the signal behavior into the reconstruction to accelerate the acquisition of MR diffusion data. JFPJ Abascal, J Parra-Robles, M Desco, (Submitted for publication) 2016.** 

We propose a novel compressed sensing method that incorporates the knowledge of the signal decay into the reconstruction (SIDER) to accelerate the acquisition of MR diffusion data by undersampling in both spatial and b-value dimensions. SIDER combines TV with a penalty function that promotes sparsity across the b-direction as follows:                              

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/SIDER_equation_1.jpg)

where Nabla is the spatial gradient, which leads to TV, F is the undersampled Fourier transform, u are the ventilation images, and M is an operator that encodes the relationship between ventilation images for consecutives values of b. This relationship can be approximated using a stretched exponential model 
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/SIDER_equation_2.jpg)

## Data 
Methods are assessed using control and COPD patient data (n=9). 

## Code
We provide MATLAB code for total variation (TV) and SIDER methods. Both methods are efficiently solved using the split Bregman implementation. 

## Summary of results ##

-**SIDER for acceleration factor x7** The following videos correspond to estimated maps of mean alevolar length for fully sampled data (top) and for SIDER method for an acceleration factor of x7 (bottom) for three controls and three COPD patients (from left to right).

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Label.jpg)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_x1_x7_SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_2_x1_x7SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_3_x1_x7SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_x1_x7SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_2_x1_x7_SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_3_x1_x7_SIDER.gif)

-**SIDER vs. zero-filling and TV for acceleration factors x1, x2, x5, x7** The following videos correspond to estimated maps of mean alevolar length for (from left to right) zero-filling, TV and SIDER methods. 

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Label2.jpg)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_3_x1x2x5x7_ZF_TV_SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_1_x1x2x5x7_z3_ZF_TV_SIDER.gif)

##  Repository files ##

The repository contains the following files:

### Data and MATLAB functions ###

- **Demo_SIDER.m:** Demo that loads data and shows how to use TV and SIDER methods 

- **SIDER.m:** SIDER reconstruction method

### Videos of results ###

#### High dose ####
- **Videos_SIDER_x7.ppt:** Videos show five slices across the volume for three control and three patient data sets


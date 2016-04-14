# compressed-sensing-diffusion-lung-MRI

This repository contains data, code and results for a novel compressed sensing method presented in the paper **Incorporation of prior knowledge of the signal behavior into the reconstruction to accelerate the acquisition of MR diffusion data. JFPJ Abascal, M Desco, J Parra-Robles, (submitted for publication) 2016.** 

The proposed method incorporates the knowledge of the signal decay into the reconstruction (SIDER) to accelerate the acquisition of MR diffusion data by undersampling in both spatial and b-value dimensions. SIDER combines total variation (TV) with a penalty function that promotes sparsity across the b-direction as follows:                              

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/SIDER_equation_1.jpg)

where Nabla is the spatial gradient, which leads to TV, F is the undersampled Fourier transform, u are the ventilation images, and M is an operator that encodes the relationship between ventilation images for consecutives values of b. This relationship can be approximated using a stretched exponential model 
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/SIDER_equation_2.jpg)

where D and alpha are estimated average values of diffusion and heterogeneity index, respectively, which can be used to estimate the mean alveolar length (Lm). The following figure shows images of ventilation for a control and a patient (top left), the signal decay (top right) and estimated maps of D, alpha and Lm (bottom).  

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Image_ControlPatient_VentilationImage_D_D_alpha_Lm.jpg)

## Data 
Methods are assessed using fully sampled diffusion datasets of three normal volunteers and three patients with COPD (n=8, two patients had two acquisitions at different sessions), available from earlier work [Parra-Robles et al., Proceedings of ISMRM 2012, 820; Parra-Robles et al., Proceedings of ISMRM 2014, 3529]. Data consisted of five slices (10 mm thick with 10 mm gap between slices), 64x64 resolution and 5 b-values (0, 1.6, 3.2, 4.8 and 6.4 s/cm2). 

## Code
We provide MATLAB code for total variation (TV) and SIDER methods. Both methods are efficiently solved using the split Bregman implementation. 

A demo file shows how to upload data, undersampled fully sampled data and to use TV and SIDER methods.   

## Summary of results ##

-**SIDER for acceleration factor x7** The following videos correspond to estimated maps of mean alveolar length for fully sampled data (top) and for SIDER method for an acceleration factor of x7 (bottom) for three controls and three COPD patients (from left to right).

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Label.jpg)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_x1_x7_SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_2_x1_x7SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_3_x1_x7SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_x1_x7SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_2_x1_x7_SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_3_x1_x7_SIDER.gif)

-**SIDER vs. zero-filling and TV for several acceleration factors** The following videos show one slice of mean alveolar length estimated at different acceleration factors (x1, x2, x4, x7) for zero-filling, TV and SIDER methods for a control (top) and COPD patient (bottom). Videos for all data sets are given in See VideosSeveralAccelerationFactors.pptx. 

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Label2.jpg)

![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/Control_1_x1x2x5x7_z3_ZF_TV_SIDER.gif)
![](https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI/blob/master/PatientCOPD_3_x1x2x5x7_ZF_TV_SIDER.gif)

##  Repository files ##

The repository contains the following files:

### Data and MATLAB functions ###

- **Demo_SIDER.m:** Demo that loads data and shows how to use TV and SIDER methods. 

- **SIDER.m:** SIDER reconstruction method.

- **TV.m:** TV reconstruction method.

- **genPDF.m:** Function to generate a pdf with polynomial variable density sampling, by Michael Lustig 2007, downloaded from (SparseMRI V0.2), http://web.stanford.edu/~mlustig/SparseMRI.html, see M. Lustig, D.L

- **genSampling_LIM.m:** Monte-carlo algorithm to generate a sampling pattern. Modified from the original function by Michael Lustig 2007

- **maxSidePeakRatio.m:** Computes the maximum sidelobe to peak ratio of point spread function for undersampled data. Used within genSampling_LIM.m

- **DataControl.mat:** Data set for a control subject

### Videos of results ###

#### Videos of all slices for fully sampled data and SIDER for x7 ####
- **Control_x1_x7_SIDER.gif** Top: Video that shows five slices across the volume of mean alveolar length for fully sampled data for a control data set. Bottom: Maps of mean alveolar length for SIDER method and for x7. 

- **Control_2_x1_x7SIDER.gif, Control_3_x1_x7SIDER.gif** Same as the previous video but for control data sets 2 and 3. 

- **PatientCOPD_x1_x7SIDER.gif** Top: Video that shows five slices across the volume of mean alveolar length for fully sampled data for a COPD patient. Bottom: Maps of mean alveolar length for SIDER method and for x7.  
 
- **PatientCOPD_2_x1_x7_SIDER.gif, PatientCOPD_3_x1_x7_SIDER.gif** Same as the previous video but for patient data sets 2 and 3. 

#### Videos of one slice for acceleration factors x1, x2, x5, x7 for all methods ####

- **Control_1_x1x2x5x7_z3_ZF_TV_SIDER.gif** From left to right, maps of mean alveolar length for ZF, TV and SIDER methods, for a control subject.

- **PatientCOPD_3_x1x2x5x7_ZF_TV_SIDER.gif** Same as previous video for a COPD patient. 

- **VideosSeveralAccelerationFactors.pptx** Same as previous video for all subjects


If you need to contact the author, please do so at jmparra@hggm.es, juanabascal78@gmail.com, desco@hggm.es

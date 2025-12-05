# GPU-enabled-Hessian-SIM
The original Hessian-SIM code has been revised to fully exploit the GPU. This work has been published in the journal, International Journal of Biomedical Imaging.

# Refer to this work
If you are using our work or are motivated by it, please use the BibTeX format below:
```plain
@article{https://doi.org/10.1155/2024/8862387,
author = {Oh, Kwangsung and Bianco, Piero R.},
title = {Facile Conversion and Optimization of Structured Illumination Image Reconstruction Code into the GPU Environment},
journal = {International Journal of Biomedical Imaging},
volume = {2024},
number = {1},
pages = {8862387},
doi = {https://doi.org/10.1155/2024/8862387},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1155/2024/8862387},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1155/2024/8862387},
abstract = {Superresolution, structured illumination microscopy (SIM) is an ideal modality for imaging live cells due to its relatively high speed and low photon-induced damage to the cells. The rate-limiting step in observing a superresolution image in SIM is often the reconstruction speed of the algorithm used to form a single image from as many as nine raw images. Reconstruction algorithms impose a significant computing burden due to an intricate workflow and a large number of often complex calculations to produce the final image. Further adding to the computing burden is that the code, even within the MATLAB environment, can be inefficiently written by microscopists who are noncomputer science researchers. In addition, they do not take into consideration the processing power of the graphics processing unit (GPU) of the computer. To address these issues, we present simple but efficient approaches to first revise MATLAB code, followed by conversion to GPU-optimized code. When combined with cost-effective, high-performance GPU-enabled computers, a 4- to 500-fold improvement in algorithm execution speed is observed as shown for the image denoising Hessian-SIM algorithm. Importantly, the improved algorithm produces images identical in quality to the original.},
year = {2024}}
```

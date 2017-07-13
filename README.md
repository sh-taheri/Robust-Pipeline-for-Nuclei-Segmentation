# Robust-Pipeline-for-Nuclei-Segmentation

Despite considerable research, nuclei segmentation is still a challenging task due noise, nonuniform illumination, and most importantly, in 2D projection images, overlapping and touching nuclei. The proposed framework is a region-based segmentation method, which consists of three major modules: i) the image is passed through a color deconvolution step to extract the desired stains; ii) then the generalized fast radial symmetry transform is applied to the image followed by non-maxima suppression to specify the initial seed points for nuclei, and their corresponding GFRS ellipses which are interpreted as the initial nuclei borders for segmentation; iii) finally, these nuclei border initial curves are evolved through the use of a statistical level-set approach along with topology preserving criteria for segmentation and separation of nuclei at the same time. The proposed method is evaluated using Hematoxylin and Eosin, and fluorescent stained images, performing qualitative and quantitative analysis, showing that the method outperforms thresholding and watershed segmentation approaches.




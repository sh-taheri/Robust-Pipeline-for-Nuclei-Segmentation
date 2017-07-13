
# Cyto-histopathological Imagery
Recently, computerized methods including automatic detection, segmentation and classification of objects in cyto-histological specimens, have drawn increased attention in the field of digital pathology as the result of developments in digital whole slide scanners and computer hardwares. Due to the essential role of nucleus in cellular functionality, automatic segmentation of cell nuclei is a fundamental prerequisite for all histological and cytological automated systems. Moreover, it substantially facilitates the segmentation of cytoplasms and the surrounding tissues.

<img src="https://github.com/sh-taheri/Robust-Pipeline-for-Nuclei-Segmentation/blob/master/images/Cyto-histopathological Imagery.png?raw=true" alt="Drawing" width="700">

# Robust-Pipeline-for-Nuclei-Segmentation

Despite considerable research, nuclei segmentation is still a challenging task due noise, nonuniform illumination, and most importantly, in 2D projection images, overlapping and touching nuclei. The proposed framework is a region-based segmentation method, which consists of three major modules: i) the image is passed through a color deconvolution step to extract the desired stains; ii) then the generalized fast radial symmetry transform is applied to the image followed by non-maxima suppression to specify the initial seed points for nuclei, and their corresponding GFRS ellipses which are interpreted as the initial nuclei borders for segmentation; iii) finally, these nuclei border initial curves are evolved through the use of a statistical level-set approach along with topology preserving criteria for segmentation and separation of nuclei at the same time. The proposed method is evaluated using Hematoxylin and Eosin, and fluorescent stained images, performing qualitative and quantitative analysis, showing that the method outperforms thresholding and watershed segmentation approaches.


# Results

<img src="https://github.com/sh-taheri/Robust-Pipeline-for-Nuclei-Segmentation/blob/master/images/Comparison1.png?raw=true" alt="Drawing" width="400"> <img src="https://github.com/sh-taheri/Robust-Pipeline-for-Nuclei-Segmentation/blob/master/images/Comparison2.png?raw=true" alt="Drawing" width="400">


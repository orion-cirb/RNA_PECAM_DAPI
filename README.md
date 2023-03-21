# RNA_PECAM_DAPI

* **Developed for:** Vianney
* **Team:** Brunet
* **Date:** March 2023
* **Software:** Fiji


### Images description

3D images taken with a x63 objective on an Airyscan

2 channels:
  1. *488:* Gene1 foci
  2. *555:* vaisseaux PECAM
  3. *405:* DAPI
  4. *642:* Gene2 foci
  
A *.roi* or *.zip* file containing ROI(s) can be provided with each image.

### Plugin description

In each ROI,
* Detect Gene1 foci with Median filtering + DoG filtering + MaxEntropy thresholding
* Detect Gene2 foci with Median filtering + DoG filtering + MaxEntropy thresholding
* Detect nuclei with cellPose (cyto2)
* Find foci inside nucleus_vessel+/-

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **CellPose** 

### Version history

Version 1 released on March 21, 2023.

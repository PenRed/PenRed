# Container description files

This folder contains the dockerfiles used to create the containerized versions of PenRed. These are the provided Dockerfiles:

* **Complete**: Creates a Fedora based image with all the required packages and software to compile and run PenRed package. 
* **Alpine**: Creates a Alpine based image with only the runtime libs and PenRed executable. Doesn't include DICOM suport.
* **Alpine-dicom**: Creates a Alpine based image with only the runtime libs and the PenRed executable. Includes DICOM suport.

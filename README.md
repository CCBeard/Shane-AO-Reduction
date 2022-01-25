# Shane-AO-Reduction


This is a public repository designed to help people who want to reduce astronomical images taken using the 3m Shane Telescope at Lick Observatory!

It is specifically designed with high contrast imaging in mind, and in particular works best when utilized with a 5 point dither sequence outlined in (cite)

The notebook "Image_reduction_plots.ipynb" is currently the main action, and most of what is needed to reduce data.

The "veto.ipynb" is required to generate a list of vetoed frames from a night of observing. That code is designed to automatically reject overexposed data (counts > 25k) or underexposed data (counts < 10k), but can be adjusted for your scientific needs

Both codes work, but this repository is still a work in progress in making it easy for other people to use. Right now, I am working on several things:

1. Reorganizing the helper functions into a .utils file that is imported, rather than having them take up space in the main file
2. Creating a detailed tutorial and a text data reduction for other people's convenience
3. Translating the notebooks into .py files in case anyone prefers that usage
4. Efficiency and annotation!

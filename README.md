# Shane-AO-Reduction



This is a public repository designed to help people who want to reduce astronomical images taken using the 3 m Shane Telescope at Lick Observatory.

It is specifically designed with high contrast imaging in mind, and in particular works best when utilized with a 5 point dither sequence outlined in [Furlan et al, 2017](https://ui.adsabs.harvard.edu/abs/2017AJ....153...71F/abstract).

The notebook "Image_reduction_plots.ipynb" is currently the main action, and most of what is needed to reduce data.

The "veto.ipynb" is required to generate a list of vetoed frames from a night of observing. That code is designed to automatically reject overexposed data (counts > 25k) or underexposed data (counts < 10k), but can be adjusted for your scientific needs.


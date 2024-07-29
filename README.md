# Streamlit-Glycosylator
[![Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://glycosylator.streamlit.app)

A Streamlit Web Application for the [Glycosylator](https://github.com/ibmm-unibe-ch/glycosylator) Python library. Using this app users can:

- atomic models for new glycans
- modify existing glycans
- glycosylate proteins and membranes
- model the shielding effect of glycans on the surface of proteins or membranes

# Running in Streamlit Cloud
The app is available on the public Streamlit Community Cloud via [glycosylator.streamlit.app](glycosylator.streamlit.app). The app's resources are limited on the cloud but should suffice for most small- to medium-size use cases. If you find the app to run out of memory
try running it locally on your machine or an HPC cluster.

## Running locally
To run the app you will need to install the following packages: 
- [Glycosylator](https://github.com/ibmm-unibe-ch/glycosylator)
- stmol
- streamlit

Then, after cloning this repository, run:
```bash
streamlit run ./src/main.py
```
Which should open a browser window with the app. 

> Note:
> The app may not display on *Safari* (streamlit has some issues with Safari compatibility)! Use Chrome or Firefox for better stability!


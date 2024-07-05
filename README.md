# Streamlit-Glycosylator

A Streamlit Web Application for the [Glycosylator](https://github.com/ibmm-unibe-ch/glycosylator) Python library. Using this app users can:

- atomic models for new glycans
- modify existing glycans
- glycosylate proteins and membranes
- model the shielding effect of glycans on the surface of proteins or membranes

## Usage
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


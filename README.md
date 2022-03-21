This repo contains notebook(s) and python code to help show how to use and access Fusion Hub data. 
FH_Hydrosat.py contains a python class which allows a user to construct an xarray dataset from a STAC search result (items), as well as download files. Please see usage in the Fusion-hub-examples.ipynb notebook, and use as needed.

To return a successful query from the STAC API, you must generate a valid token and store it in stac-token.txt
A token can be generated by logging in to the token server at https://fusion-stac.hydrosat.com/token

This notebook and additional documentation can be found at http://fusion.hydrosat.com/docs

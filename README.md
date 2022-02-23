This repo contains notebook(s) and python code to help show how to use and access Fusion Hub data. 
FH_Hydrosat.py contains a python class which allows a user to construct an xarray dataset from a STAC search result (items), as well as download files.

To run, you must generate a valid token and store it in stac-token.txt
I've chosen this in lieu of storing credentials in the notebook as the token will expire and eventually we will have an endpoint to generate one with a GET request.
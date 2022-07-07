#!/usr/bin/env python

import requests as req

# web-site's file processing is done onload().  This script
# accesses the site, so any .gz file uploaded from a workflows-
# run is processed before a new file is uploaded.

try:
    resp = req.request(method='GET', url="http://plot.joelatessa.com")
except Exception:
    print('could not access site... continuing')
    pass

# show in terminal for testing
# print(resp.text)

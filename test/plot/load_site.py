#!/usr/bin/env python

import requests as req

resp = req.request(method='GET', url="http://plot.joelatessa.com")

# show in terminal for testing
# print(resp.text)

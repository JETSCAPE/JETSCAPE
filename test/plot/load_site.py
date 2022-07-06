#!/usr/bin/env python

import requests as req

resp = req.request(method='GET', url="http://plot.joelatessa.com")
print(resp.text)

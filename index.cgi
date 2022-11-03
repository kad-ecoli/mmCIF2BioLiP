#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
from datetime import datetime
import subprocess

rootdir=os.path.dirname(os.path.abspath(__file__))

html_header=""
html_footer=""
if os.path.isfile(rootdir+"/index.html"):
    fp=open(rootdir+"/index.html")
    txt=fp.read()
    fp.close()
    html_header=txt.split('<!-- CONTENT REFRESH -->')[0]
    html_footer=txt.split('<!-- CONTENT END -->')[-1]

print("Content-type: text/html\n")
if len(html_header):
   print(html_header)
else:
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/BioLiP1.png ></br>
<p><a href=.>[Back to Home]</a></p>
''')

filename=rootdir+"/data/index.txt"
if os.path.isfile(filename):
    fp=open(filename,'r')
    print(fp.read())
    fp.close()

if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

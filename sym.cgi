#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip

rootdir=os.path.dirname(os.path.abspath(__file__))

form = cgi.FieldStorage()
print("Content-type: text/html\n")
print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP:Ligand Information</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/BioLiP1.png ></br>
<p><a href=.>[Back to Home]</a></p>
<style>
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
}
</style>
<table>
''')

code=form.getfirst("code",'').upper()
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
txt_table=''
for line in fp.read().splitlines():
    items=line.split('\t')
    if line.startswith('#'):
        items[0]=items[0][1:]
        txt_table+="<tr><th>"+"</th><th>".join(items)+"</th></tr>"
    if code!=items[0]:
        continue
    for i in range(len(items)):
        if "; " in items[i]:
            items[i]=items[i].replace("; ",";<br>")
    txt_table+="<tr><td>"+"</td><td>".join(items)+"</td></tr>"
    res=items[0]
    svg="https://cdn.rcsb.org/images/ccd/labeled/%s/%s.svg"%(res[0],res)
    print("<p><ul><a href=%s target=_blank><img src=%s width=300></a><br>"%(svg,svg))
    print("View <a href=https://rcsb.org/ligand/%s target=_blank>%s</a> at the PDB database</ul></p>"%(
        res,res))
print(txt_table) 
print('''
</table>
</body>
</html>
''')

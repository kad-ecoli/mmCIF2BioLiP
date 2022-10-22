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
header_list=[]
for line in fp.read().splitlines():
    items=line.split('\t')
    if line.startswith('#'):
        items+=["image"]
        print("<tr><th>"+"</th><th>".join(items)+"</th></tr>")
    if code!=items[0]:
        continue
    for i in range(len(items)):
        if "; " in items[i]:
            items[i]=items[i].replace("; ",";<br>")
    res=items[0]
    svg="https://cdn.rcsb.org/images/ccd/labeled/%s/%s.svg"%(res[0],res)
    items.append(
        "<a href=%s target=_blank><img src=%s width=150>"%(svg,svg))
    items[0]="<a href=https://www.rcsb.org/ligand/%s target=_blank>%s</a>"%(
        res,res)
    print("<tr><td>"+"</td><td>".join(items)+"</td></tr>")
            



print('''
</table>
</body>
</html>
''')

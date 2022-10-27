#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip

rootdir=os.path.dirname(os.path.abspath(__file__))

html_header=""
html_footer=""
if os.path.isfile(os.path.join(rootdir,"index.html")):
    fp=open(os.path.join(rootdir,"index.html"))
    txt=fp.read()
    fp.close()
    html_header=txt.split('<!-- CONTENT START -->')[0]
    html_footer=txt.split('<!-- CONTENT END -->')[-1]

form = cgi.FieldStorage()
print("Content-type: text/html\n")
if len(html_header):
   print(html_header)
else:
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP:Ligand Information</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/BioLiP1.png ></br>
<p><a href=.>[Back to Home]</a></p>
''')

print('''
<style>
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
}
</style>
<table>
''')

code=form.getfirst("code",'')
if code in ["peptide","rna","dna"]:
    print('''
<tr><th>New BioLiP code for polymer</th><th>Old BioLiP code for Polymer</th></tr>
<tr><td>peptide</td><td>III</td></tr>
<tr><td>rna</td><td>NUC</td></tr>
<tr><td>dna</td><td>NUC</td></tr>
''')
elif code:
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
        print("<p><ul><a href=%s target=_blank><img src=%s alt='' width=300></a><br>"%(svg,svg))
        print("View <a href=https://rcsb.org/ligand/%s target=_blank>%s</a> at the PDB database</ul></p>"%(
        res,res))
    print(txt_table)
else:
    fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
    ligand_list=[]
    for line in fp.read().splitlines()[1:]:
        ligand_list.append(line.split('\t')[0])
    fp.close()
    import random
    code=random.choice(ligand_list)
    
    print("No ligand code provided. You may <a href=%s?code=%s>[browse a random ligand]</a>"%(
        os.path.basename(__file__),code
    ))
print("</table>")

if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip

rootdir=os.path.dirname(os.path.abspath(__file__))

html_header=""
html_footer=""
if os.path.isfile(rootdir+"/index.html"):
    fp=open(rootdir+"/index.html")
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
<title>BioLiP</title>
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
            elif i==2:
                inchi_list=[]
                for j in range(0,len(items[2]),20):
                    start=j
                    end  =j+20
                    inchi_list.append(items[2][start:end])
                items[2]='<br>'.join(inchi_list)
            elif i==4:
                smiles_list=[]
                for j in range(0,len(items[4]),20):
                    start=j
                    end  =j+20
                    smiles_list.append(items[4][start:end])
                items[4]='<br>'.join(smiles_list)
                items[4].replace('; ',';<br>')
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

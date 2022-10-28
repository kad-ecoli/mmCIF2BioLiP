#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import subprocess

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
page=form.getfirst("page",'')
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

# pdb recCha => resolution ec go uniprot pubmed
cmd="zcat %s/data/pdb_all.tsv.gz|cut -f1-3,6-9"%rootdir
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
stdout=stdout.decode()
chain_dict=dict()
for line in stdout.splitlines()[1:]:
    items=line.split('\t')
    key='\t'.join(items[:2])
    value='\t'.join(items[2:])
    chain_dict[key]=value

ligand_dict=dict()
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    ccd  =items[0]
    name =items[1]
    ligand_dict[ccd]=name
fp.close()


fp=gzip.open(rootdir+"/data/lig_all.tsv.gz",'rt')
lines=fp.read().splitlines()[1:]
fp.close()

totalNum=len(lines)

print('''
Download all results in tab-seperated text for 
<a href=data/pdb_all.tsv.gz>all %d receptors</a> and
<a href=data/lig_all.tsv.gz>all %d receptor-ligand interactions</a>.<br>
Click <strong>Site #</strong> to view the binding site.<br>
Hover over <strong>Ligand</strong> to view the full ligand name.<br>
Hover over <strong>GO terms</strong> to view all GO terms.
<p></p>
'''%(len(chain_dict),totalNum))

pageLimit=200
totalPage=1+int(totalNum/pageLimit)
if not page:
    page=1
elif page=="last":
    page=totalPage
else:
    page=int(page)
if page<1:
    page=1
elif page>totalPage:
    page=totalPage

print('''<center> 
<a class="hover" href="?&page=1">&lt&lt</a>
<a class="hover" href="?&page=%d">&lt</a>
'''%(page-1))
for p in range(page-10,page+11):
    if p<1 or p>totalPage:
        continue
    elif p==page:
        print(' %d '%(p))
    else:
        print('<a class="hover" href="?&page=%d">%d</a>'%(p,p))
print('''
<a class="hover" href="?&page=%d">&gt</a>
<a class="hover" href="?&page=last">&gt&gt</a>
<form name="pform" action="alllig.cgi">Go to page <select name="page" onchange="this.form.submit()">
'''%(page+1))
for p in range(1,totalPage+1):
    if p==page:
        print('<option value="%d" selected="selected">%d</option>'%(p,p))
    else:
        print('<option value="%d">%d</option>'%(p,p))
print("</select></form></center><br>")
     
print('''  
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </srong></th>
    <th width=10% ALIGN=center><strong> PDB<br>(Resolution &#8491;) </srong></th>
    <th width=10% ALIGN=center><strong> Site # </srong></th>
    <th width=10% ALIGN=center><strong> Ligand </srong> </th>           
    <th width=10% ALIGN=center><strong> EC number </srong> </th>           
    <th width=10% ALIGN=center><strong> GO terms </srong> </th>           
    <th width=15% ALIGN=center><strong> UniProt </srong> </th>           
    <th width=10% ALIGN=center><strong> PubMed </srong> </th>           
    <th width=20% ALIGN=left> <strong> Binding affinity</srong> </th>           
</tr><tr ALIGN=center>
''')

sort_line=[]
for l,line in enumerate(lines):
    items     =line.split('\t')
    pdb       =items[0]
    recCha    =items[1]
    bs        =items[2]
    ccd       =items[3]
    name      =""
    if ccd in ligand_dict:
        name=';<br>'.join(ligand_dict[ccd].split(';'))
    ligCha    =items[4]
    manual    =items[5]
    moad      =items[6]
    pdbbind   =items[7]
    bindingdb=items[8]

    resolution="N/A"
    ec        ="N/A"
    go        ="N/A"
    uniprot   ="N/A"
    pubmed    ="N/A"
    key       =pdb+'\t'+recCha
    if key in chain_dict:
        items     =chain_dict[key].split('\t')
        resolution=items[0]
        ec        =items[1]
        go        =items[2]
        uniprot   =items[3]
        pubmed    =items[4]
    
   
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><a href="qsearch.cgi?lig3=%s" target="_blank">%s</td>
    <td>%s</td>
    <td ALIGN=left><a href="sym.cgi?code=%s" target="_blank">%s</td>
</tr>
'''%(bgcolor,l+1,ccd,ccd,'',ccd,name))

for l in range(pageLimit*(page-1),pageLimit*page+1):
    if l>=totalNum:
        continue


print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

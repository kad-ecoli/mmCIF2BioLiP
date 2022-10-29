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
page =form.getfirst("page",'')
lig3 =form.getfirst("lig3",'').lower()
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

# pdb recCha => resolution (csa csa_renumbered) ec go uniprot pubmed
fp=gzip.open(rootdir+"/data/pdb_all.tsv.gz",'rt')
chain_dict=dict()
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    chain=':'.join(items[:2])
    chain_dict[chain]=items[2:3]+items[5:]
fp.close()

if not lig3 in ["rna","dna","peptide"]:
    print('''
<li><a href=polymer.cgi?lig3=rna>Browse RNA ligands</a></li>
<li><a href=polymer.cgi?lig3=dna>Browse DNA ligands</a></li>
<li><a href=polymer.cgi?lig3=peptide>Browse peptide ligands</a></li>
    ''')
    if len(html_footer):
        print(html_footer)
    else:
        print("</body> </html>")
    exit()

fasta_dict=dict()
fp=gzip.open("%s/data/%s.fasta.gz"%(rootdir,lig3),'rt')
chain_list=[]
for block in fp.read().split('>')[1:]:
    header,sequence=block.splitlines()
    chain=header.split()[0]
    fasta_dict[chain]=sequence
    chain_list.append(chain)
fp.close()
chain_list.sort()

clust_dict=dict()
fp=gzip.open("%s/data/%s_nr.fasta.clust.gz"%(rootdir,lig3),'rt')
for line in fp.read().splitlines():
    rep,mem=line.split('\t')
    mem_list=[rep]+mem.split(',')
    for mem in mem_list:
        clust_dict[mem]=[m.replace('_%s_'%lig3,':') for m in mem_list if m!=mem]
fp.close()

cmd="zcat %s/data/lig_all.tsv.gz|grep -P '^\w+\\t\w+\\tBS\d+\\t%s\\t'"%(rootdir,lig3)
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
lig2rec_dict=dict()
for line in stdout.decode().splitlines():
    items=line.split('\t')
    pdb  =items[0]
    ligCha=items[4]
    ligKey='_'.join((pdb,lig3,ligCha))
    if not ligKey in lig2rec_dict:
        lig2rec_dict[ligKey]=[]
    lig2rec_dict[ligKey].append(items)

totalNum=len(fasta_dict)

print('''
Download all <a href=data/%s.fasta.gz>%d %s sequences</a>.<br>
Resolution -1.00 means the resolution is unavailable, e.g., for NMR structures.
Click <strong>Receptor chain</strong> to view the binding site structure.
Hover over <strong>Receptor chain</strong> to view the binding residues.
The sequence is converted from residues with experimentally determined coordinates in the structure; residues not observed in the 3D structure are excluded.
<p></p>
'''%(lig3,totalNum,lig3  if lig3=="peptide" else lig3.upper()))

pageLimit=100
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
<a class="hover" href="?page=1&lig3=%s">&lt&lt</a>
<a class="hover" href="?page=%d&lig3=%s">&lt</a>
'''%(lig3,page-1,lig3))
for p in range(page-10,page+11):
    if p<1 or p>totalPage:
        continue
    elif p==page:
        print(' %d '%(p))
    else:
        print('<a class="hover" href="?page=%d&lig3=%s">%d</a>'%(p,lig3,p))
print('''
<a class="hover" href="?&page=%d&lig3=%s">&gt</a>
<a class="hover" href="?&page=last&lig3=%s">&gt&gt</a>
<form name="pform" action="polymer.cgi">Go to page
<input type="hidden" name="lig3" value="%s" size=0>
<select name="page" onchange="this.form.submit()">
'''%(page+1,lig3,lig3,lig3))
for p in range(1,totalPage+1):
    if p==page:
        print('<option value="%d" selected="selected">%d</option>'%(p,p))
    else:
        print('<option value="%d">%d</option>'%(p,p))
print("</select></form></center><br>")
     
print('''  
<table border="0" align=center width=100% style="table-layout: fixed;">    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=10% ALIGN=center><strong> Ligand chain<br>(Resolution &#8491;) </strong></th>
    <th width=10% ALIGN=center><strong> Receptor chain<br>(Site #) UniProt</strong></th>
    <th width=65% ALIGN=center><strong>
''')
print(lig3.upper() if lig3!="peptide" else "Peptide")
print('''
sequence </strong> </th>           
    <th width=10% ALIGN=center><strong> PubMed</strong></th>
</tr><tr ALIGN=center>
''')

for l in range(pageLimit*(page-1),pageLimit*page):
    if l>=totalNum:
        continue
    ligKey=chain_list[l]
    items =lig2rec_dict[ligKey][0]
    sequence=fasta_dict[ligKey]
    if ligKey in clust_dict and len(clust_dict[ligKey]):
        sequence="(Identical to "+'; '.join(clust_dict[ligKey]
            )+')<br>'+sequence
    pdb,ligCha=ligKey.split('_'+lig3+'_')

    bs_list=[]
    reso      ="N/A"
    pubmed    ="N/A"
    for i in range(len(lig2rec_dict[ligKey])):
        items     =lig2rec_dict[ligKey][i]
        recCha    =items[1]
        bs        =items[2]
        resOrig   =items[6]
        resRenu   =items[7]
        key       =pdb+':'+recCha
        uniprot   ="N/A"
        if key in chain_dict:
            items     =chain_dict[key]
            reso      =items[0]
            uniprot   =items[3]
            pubmed    =items[4]
            if uniprot:
                uniprot=','.join(["<a href=https://uniprot.org/uniprot/%s target=_blank>%s</a>"%(u,u) for u in uniprot.split(',')])
            else:
                uniprot="N/A"
            if pubmed:
                pubmed="<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(pubmed,pubmed)
            else:
                pubmed="N/A"

        bs_list.append('''
<span title="%s"><a href="getaid.cgi?pdb=%s&chain=%s&bs=%s" target=_blank>%s (%s)</a></span> %s
'''%(resOrig,pdb,recCha,bs,recCha,bs,uniprot))
    
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><a href="pdb.cgi?pdb=%s&chain=%s" target=_blank>%s:%s</a> (%s)</td>
    <td>%s</td>
    <td style="word-wrap: break-word">%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    l+1,
    pdb,ligCha,pdb,ligCha,reso,
    '<br>'.join(bs_list),
    sequence,
    pubmed,
    ))


print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

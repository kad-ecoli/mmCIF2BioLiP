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

#### read cgi parameters ####

form = cgi.FieldStorage()
page =form.getfirst("page",'')
order=form.getfirst("order",'').lower()
if not order:
    order="pdbid"
#### read database data ####
# pdb:recCha => resolution ec go uniprot pubmed
fp=gzip.open(rootdir+"/data/pdb_all.tsv.gz",'rt')
chain_dict=dict()
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    chain_dict[':'.join(items[:2])]=items[2:3]+items[5:]
fp.close()

ligand_dict=dict()
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    ccd  =items[0]
    name =items[-1]
    ligand_dict[ccd]=name
fp.close()


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

enzyme_dict=dict()
fp=gzip.open(rootdir+"/data/enzyme.tsv.gz",'rt')
for line in fp.read().splitlines():
    e,name=line.split('\t')
    enzyme_dict[e]=name
fp.close()

go2name_dict=dict()
fp=gzip.open(rootdir+"/data/go2name.tsv.gz",'rt')
for line in fp.read().splitlines():
    g,a,name=line.split('\t')
    go2name_dict[g]='('+a+') '+name
fp.close()

#### parse page ####
pageLimit=200


fp=gzip.open(rootdir+"/data/lig_all.tsv.gz",'rt')
lines=fp.read().splitlines()[1:]
fp.close()

totalNum=len(lines)

print('''
Download all results in tab-seperated text for 
<a href=data/pdb_all.tsv.gz>%d receptors</a> and
<a href=data/lig_all.tsv.gz>%d receptor-ligand interactions</a>.<br>
<li>Click <strong>PDB</strong> to view the structure at the RCSB PDB database.
Resolution -1.00 means the resolution is unavailable, e.g., for NMR structures.</li>
<li>Click <strong>Site #</strong> to view the binding site structure.
Hover over <strong>Site #</strong> to view the binding residues.</li>
<li>Hover over <strong>Ligand</strong> to view the full ligand name.</li>
<li>Hover over <strong>EC number</strong> to view the full name of enzymatic activity.</li>
<li>Hover over <strong>GO terms</strong> to view all GO terms.
Click <strong>GO terms</strong> to view the GO annotations for the UniProt protein associated with the PDB chain</li>
<p></p>
'''%(len(chain_dict),totalNum))

print('''
<form name="sform" action="browse.cgi">
Sort results by
<select name="order" onchange="this.form.submit()">
    <option value="pdbid">PDB ID</option>
    <option value="lig3">Ligand ID</option>
    <option value="uniprot">UniProt ID</option>
    <option value="reso">Resolution</option>
</select>
</form>
'''.replace('value="%s"'%order,
            'value="%s" selected="selected"'%order))

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
<form name="pform" action="browse.cgi">Go to page <select name="page" onchange="this.form.submit()">
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
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=12% ALIGN=center><strong> PDB<br>(Resolution &#8491;) </strong></th>
    <th width=5%  ALIGN=center><strong> Site # </strong></th>
    <th width=10% ALIGN=center><strong> Ligand </strong> </th>           
    <th width=10% ALIGN=center><strong> EC number </strong> </th>           
    <th width=16% ALIGN=center><strong> GO terms </strong> </th>           
    <th width=10% ALIGN=center><strong> UniProt </strong> </th>           
    <th width=10% ALIGN=center><strong> PubMed </strong> </th>           
    <th width=22% ALIGN=center><strong> Binding<br>affinity</strong> </th>           
</tr><tr ALIGN=center>
''')

sort_line=[]
for l,line in enumerate(lines):
    items=line.split('\t')
    chain=':'.join(items[:2])
    if order=="lig3":
        sort_line.append((items[3],items))
    elif order=="uniprot":
        if not chain in chain_dict:
            continue
        uniprot=chain_dict[chain][3]
        sort_line.append((uniprot,items))
    elif order=="reso":
        if not chain in chain_dict:
            continue
        reso=chain_dict[chain][0]
        sort_line.append((reso,items))
    elif order=="pdbid":
        sort_line.append((chain,items))
sort_line.sort()

for l in range(pageLimit*(page-1),pageLimit*page):
    if l>=totalNum:
        continue
    items=sort_line[l][1]
    pdb       =items[0]
    recCha    =items[1]
    bs        =items[2]
    ccd       =items[3]
    name      =""
    if ccd in ligand_dict:
        name=';\n'.join(ligand_dict[ccd].split(';'))
    ligCha    =items[4]
    ligIdx    =items[5]
    resOrig   =items[6]
    resRenu   =items[7]
    manual    =items[8]
    moad      =items[9]
    pdbbind   =items[10]
    bindingdb =items[11]
    
    affinity  =''
    if manual:
        affinity="Manual survey: "+manual.replace(',',', ')+"<br>"
    if moad:
        affinity+="<a href=http://bindingmoad.org/pdbrecords/index/%s target=_blank>MOAD</a>: %s<br>"%(pdb,moad.replace(',',', '))
    if pdbbind:
        affinity+='<a href="http://pdbbind.org.cn/quickpdb.php?quickpdb=%s" target=_blank>PDBbind:</a> %s<br>'%(pdb,pdbbind.replace(',',', '))
    if bindingdb:
        affinity+="BindingDB: "+bindingdb.replace(',',', ')+"<br>"
    if affinity:
        affinity=affinity[:-4]

    key       =pdb+':'+recCha

    reso      ="N/A"
    ec        ="N/A"
    go        ="N/A"
    uniprot   ="N/A"
    pubmed    ="N/A"
    if key in chain_dict:
        items     =chain_dict[key]
        reso      =items[0]
        ec        =items[1]
        go        =items[2]
        uniprot   =items[3]
        pubmed    =items[4]
        if ec:
            ec_list=ec.split(',')
            ec=''
            for e in ec_list:
                if ec:
                    ec+='<br>'
                if not e in enzyme_dict:
                    ec+="<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(e,e)
                else:
                    ec+='<a href=https://enzyme.expasy.org/EC/%s target=_blank><span title="%s">%s</span></a>'%(e,enzyme_dict[e],e)
        else:
            ec="N/A"
        if go:
            go_list=["GO:"+g for g in go.split(',')]

            go='<span title="'
            for g in go_list:
                if g in go2name_dict:
                    go+=g+' '+go2name_dict[g]+'\n'
                else:
                    go+=g+'\n'
            go=go[:-1]+'">'+go_list[0]+" ...</span>"
            if uniprot:
                go='<a href="https://ebi.ac.uk/QuickGO/annotations?geneProductId=%s" target=_blank>%s</a>'%(uniprot.split(',')[0],go)
        else:
            go="N/A"
        if uniprot:
            uniprot=','.join(["<a href=https://uniprot.org/uniprot/%s target=_blank>%s</a>"%(u,u) for u in uniprot.split(',')])
        else:
            uniprot="N/A"
        if pubmed:
            pubmed="<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(pubmed,pubmed)
        else:
            pubmed="N/A"
   
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    ccd_http=ccd
    if name:
        ccd_http='<span title="%s">%s</span>'%(name,ccd)
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><a href="https://rcsb.org/structure/%s" target=_blank>%s:%s</a> (%s)</td>
    <td><span title="%s"><a href="pdb.cgi?pdb=%s&chain=%s&bs=%s" target=_blank>%s</span></td>
    <td><a href="sym.cgi?code=%s" target=_blank>%s</a></td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    l+1,
    pdb,pdb,recCha,reso,
    resOrig,pdb,recCha,bs,bs,
    ccd,ccd_http,
    ec,
    go,
    uniprot,
    pubmed,
    affinity,
    ))


print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

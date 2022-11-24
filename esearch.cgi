#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import textwrap

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
page =form.getfirst("page",'').strip().strip("'")
if not page:
    page='1'
elif page=='0':
    page='1'
order=form.getfirst("order",'').lower().strip().strip("'")
if not order:
    order="pdbid"
pdbid=form.getfirst("pdbid",'').lower().strip().strip("'")
chain=form.getfirst("chain",'').strip().strip("'")
uniprot=form.getfirst("uniprot",'').upper().strip().strip("'")
ecn    =form.getfirst("ecn",'').strip().strip("'")
got    =form.getfirst("got",'').upper().strip().strip("'")
if got:
    got=got.split()[0]
    if got.startswith("GO:"):
        got=got[3:]
pubmed =form.getfirst("pubmed",'').strip("'")
outfmt =form.getfirst("outfmt",'').strip().strip("'")

para_list=[]
if order:
    para_list.append("order=%s"%order)
if pdbid:
    para_list.append("pdbid=%s"%pdbid)
if chain:
    para_list.append("chain=%s"%chain)
if uniprot:
    para_list.append("uniprot=%s"%uniprot)
if ecn:
    para_list.append("ecn=%s"%ecn)
if got:
    para_list.append("got=%s"%got)
if pubmed:
    para_list.append("pubmed=%s"%pubmed)
para='&'.join(para_list)

#### read database data ####
if outfmt!='txt':
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

got2chain_set=set()
if got:
    got2chain_list=[]
    fp=gzip.open(rootdir+"/data/pdb_go.tsv.gz",'rt')
    for line in fp.read().splitlines():
        pdb,recCha,go_line=line.split('\t')
        if got in go_line.split(','):
            got2chain_list.append(pdb+':'+recCha)
    fp.close()
    got2chain_set=set(got2chain_list)

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

pdb2name_dict=dict()
fp=gzip.open(rootdir+"/data/title.tsv.gz",'rt')
for line in fp.read().splitlines():
    p,name=line.split('\t')[:2]
    pdb2name_dict[p]=name
fp.close()

sprot_dict=dict()
fp=gzip.open(rootdir+"/data/uniprot_sprot.tsv.gz",'rt')
for line in fp.read().splitlines():
    u,name,gn=line.split('\t')
    if gn:
        name+=" (Gene Name="+gn+")"
    sprot_dict[u]=name
fp.close()

#### parse page ####
pageLimit=200
html_txt=''
sort_line=[]
# pdb:recCha => resolution csa csa_renumbered ec go uniprot pubmed
fp=gzip.open(rootdir+"/data/ec_all.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
#for line in stdout.decode().splitlines():
    items=line.split('\t')
    pdb       =items[0]
    recCha    =items[1]
    reso      =items[2]
    csaOrig   =items[3]
    csaRenu   =items[4]
    ec        =items[5]
    go        =items[6]
    accession =items[7]
    pmid      =items[8]
    if pdbid and pdb!=pdbid:
        continue
    if chain and recCha!=chain:
        continue
    if got and not pdb+':'+recCha in got2chain_set:
        continue
    if uniprot and not uniprot in accession.split(','):
        continue
    if ecn and ecn!='0' and not ecn in ec.split(','):
        continue
        go=items[-3]
    if pubmed and pmid!=pubmed:
        continue

    items=(pdb,recCha,reso,csaOrig,csaRenu,ec,go,accession,pmid)
    if outfmt=='txt':
        html_txt+='\t'.join(items)+'\n'
    else:
        if order=="reso":
            sort_line.append((reso,items))
        elif order=="ecn":
            sort_line.append((ecn,items))
        elif order=="uniprot":
            sort_line.append((accession,items))
        else:
            sort_line.append((pdb+recCha,items))

if outfmt=="txt":
    print("Content-type: text/plain\n")
    print(html_txt)
    exit()

sort_line.sort()
totalNum=len(sort_line)
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

for l in range(totalNum):
    if l<pageLimit*(int(page)-1) or l>=pageLimit*(int(page)):
        continue
    items    =sort_line[l][1]
    
    pdb       =items[0]
    recCha    =items[1]
    reso      =items[2]
    csaOrig   =items[3]
    csaRenu   =items[4]
    ec        =items[5]
    go        =items[6]
    accession =items[7]
    pmid      =items[8]
    
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
        if accession:
            go='<a href="https://ebi.ac.uk/QuickGO/annotations?geneProductId=%s" target=_blank>%s</a>'%(accession.split(',')[0],go)
    else:
        go="N/A"
    if accession:
        #accession='<br>'.join(["<a href=https://uniprot.org/uniprot/%s target=_blank>%s</a>"%(a,a) for a in accession.split(',')])
        accession_list=[]
        for a in accession.split(','):
            name=''
            if a in sprot_dict:
                name=sprot_dict[a].replace('"','')
            a="<a href=https://uniprot.org/uniprot/%s target=_blank>%s</a>"%(a,a)
            if name:
                a='<span title="%s">%s</span>'%(name,a)
            accession_list.append(a)
        accession='<br>'.join(accession_list)
    else:
        accession="N/A"
    if pmid:
        pmid="<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(pmid,pmid)
    else:
        pmid="N/A"
    
    csa=''
    if csaRenu:
        csa='''<span title="%s"><a href="https://www.ebi.ac.uk/thornton-srv/m-csa/search/?s=%s" target=_blank>%s</a></span>'''%(
            csaRenu, pdb,csaOrig)

    reso="("+reso+")"
    title=''
    if pdb in pdb2name_dict:
        title=pdb2name_dict[pdb]
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    html_txt+='''
<tr %s ALIGN=center>
    <td><a href="pdb.cgi?pdb=%s&chain=%s" target=_blank>%d</a></td>
    <td><span title="%s"><a href="https://rcsb.org/structure/%s" target=_blank>%s:%s</a> %s</span></td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    pdb,recCha,l+1,
    title,pdb,pdb,recCha,reso,
    csa,
    ec,
    go,
    accession,
    pmid,
    )
fp.close()

print('''
<h4>Search Enzyme</h4>
<FORM id="form1" name="form1" METHOD="POST" ENCTYPE="MULTIPART/FORM-DATA" ACTION="esearch.cgi">
    PDB ID
    <span title="Example: 1a69">
    <input type="text" name="pdbid" value="" placeholder="%s" size=9>
    </span>
    &nbsp;&nbsp;
    
    Chain ID
    <span title="Example: A">
    <input type="text" name="chain" value="" placeholder="%s" size=9>
    </span>
    &nbsp;&nbsp;
    
    UniProt accession
    <span title="Example: P0ABP8">
    <input type="text" name="uniprot" value="" placeholder="%s" size=9>
    </span>
    &nbsp;&nbsp;
    
    EC number
    <span title="Example: 2.4.2.1">
    <input type="text" name="ecn" value="" placeholder="%s" size=9>
    </span>
    &nbsp;&nbsp;
    
    GO term
    <span title="Example: 0004731, or GO:0004731">
    <input type="text" name="got" value="" placeholder="%s" size=9>
    </span>
    &nbsp;&nbsp;
   
    PubMed ID
    <span title="Example: 9653038">
    <input type="text" name="pubmed" value="" placeholder="%s" size=9>
    </span>
    &nbsp;&nbsp;

    <INPUT TYPE="submit" VALUE="Submit">
    &nbsp;
    <INPUT TYPE="reset" VALUE="Clear">
</FORM>
<p></p>
'''%(pdbid,chain,uniprot,ecn,got,pubmed))

print('''
<h4>Browse Enzyme</h4>
This page include enzymes with and without ligand interactions. Download all results in tab-seperated text for 
<a href="?outfmt=txt&%s" download="BioLiP.txt">%d</a> enzymes, whose format is explained at <a href="download/readme_enzyme.txt">readme_enzyme.txt</a>.<br>
<li>Click <strong>#</strong> to view the structure.</li>
<li>Hover over <strong>PDB</strong> to view the title of the structure.
Click <strong>PDB</strong> to view the structure at the RCSB PDB database.
Resolution -1.00 means the resolution is unavailable, e.g., for NMR structures.</li>
'''%(para,totalNum))
print('''
<li>Hover over <strong>EC number</strong> to view the full name of enzymatic activity.</li>
<li>Hover over <strong>GO terms</strong> to view all GO terms.
<li>Hover over <strong>UniProt</strong> to view the protein name.</li>
''')


print(('''<p></p>
<form name="sform" action="esearch.cgi">
Sort results by
<select name="order" onchange="this.form.submit()">
    <option value="pdbid">PDB ID</option>
    <option value="ecn">EC number</option>
    <option value="uniprot">UniProt ID</option>
    <option value="reso">Resolution</option>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=uniprot value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=pubmed  value='%s'>
</form>'''%(pdbid,uniprot,ecn,got,pubmed)
).replace('value="%s"'%order,
          'value="%s" selected="selected"'%order))


print('''<center> 
<a class='hover' href='?&page=1&%s'>&lt&lt</a>
<a class='hover' href='?&page=%d&%s'>&lt</a>
'''%(para,page-1,para))
for p in range(page-10,page+11):
    if p<1 or p>totalPage:
        continue
    elif p==page:
        print(' %d '%(p))
    else:
        print('''<a class='hover' href='?&page=%d&%s'>%d</a>'''%(p,para,p))
print('''
<a class='hover' href='?&page=%d&%s'>&gt</a>
<a class='hover' href='?&page=last&%s'>&gt&gt</a>
<form name="pform" action="ssearch.cgi">Go to page <select name="page" onchange="this.form.submit()">
'''%(page+1,para,para))
for p in range(1,totalPage+1):
    if p==page:
        print('<option value="%d" selected="selected">%d</option>'%(p,p))
    else:
        print('<option value="%d">%d</option>'%(p,p))
print('''</select>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=uniprot value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=pubmed  value='%s'>
</form></center><br>'''%(pdbid,uniprot,ecn,got,pubmed))

print('''  
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=7% ALIGN=center><strong> # </strong></th>
    <th width=14% ALIGN=center><strong> PDB<br>(Resolution &#8491;) </strong></th>
    <th width=31%  ALIGN=center><strong> Catalytic<br>Site </strong></th>
    <th width=10% ALIGN=center><strong> EC number </strong> </th>           
    <th width=18% ALIGN=center><strong> GO terms </strong> </th>           
    <th width=10% ALIGN=center><strong> UniProt </strong> </th>           
    <th width=10% ALIGN=center><strong> PubMed </strong> </th>           
</tr><tr ALIGN=center>
''')
print(html_txt)
print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

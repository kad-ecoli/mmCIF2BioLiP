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
page =form.getfirst("page",'').strip().strip("'")
if not page:
    page='1'
elif page=='0':
    page='1'
pdbid=form.getfirst("pdbid",'').lower().strip().strip("'")
chain=form.getfirst("chain",'').strip().strip("'")
lig3 =form.getfirst("lig3",'').strip().strip("'")
if not lig3:
    lig3=form.getfirst("code",'').strip().strip("'")
uniprot=form.getfirst("uniprot",'').upper().strip().strip("'")
ecn    =form.getfirst("ecn",'').strip().strip("'")
got    =form.getfirst("got",'').upper().strip().strip("'")
if got:
    got=got.split()[0]
    if got.startswith("GO:"):
        got=got[3:]
ligname=form.getfirst("ligname",'').upper().strip().strip("'")
pubmed =form.getfirst("pubmed",'').strip("'")
baff   =form.getfirst("baff",'').strip().strip("'")
outfmt =form.getfirst("outfmt",'').strip().strip("'")

para_list=[]
if pdbid:
    para_list.append("pdbid='%s'"%pdbid)
if chain:
    para_list.append("chain='%s'"%chain)
if lig3:
    para_list.append("lig3='%s'"%lig3)
elif ligname:
    para_list.append("ligname='%s'"%ligname)
if uniprot:
    para_list.append("uniprot='%s'"%uniprot)
if ecn:
    para_list.append("ecn='%s'"%ecn)
if got:
    para_list.append("got='%s'"%got)
if baff:
    para_list.append("baff='baff'")
if pubmed:
    para_list.append("pubmed='%s'"%pubmed)
para='&'.join(para_list)

#### read database data ####
# pdb:recCha => resolution csa csa_renumbered ec go uniprot pubmed
fp=gzip.open(rootdir+"/data/pdb_all.tsv.gz",'rt')
chain_dict=dict()
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    chain_dict[':'.join(items[:2])]=items[2:]
fp.close()

ligand_dict=dict()
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    ccd  =items[0]
    name =items[-1]
    ligand_dict[ccd]=name
fp.close()
lig_set=set()
if lig3:
    lig_set=set([lig3])
    if lig3=="metal" or lig3=="regular":
        fp=gzip.open(rootdir+"/data/metal.tsv.gz",'rt')
        metal_set=set([line.split()[0] for line in fp.read().splitlines()])
        if lig3=="metal":
            lig_set=metal_set
        else:
            lig_set=set([ccd for ccd in ligand_dict if not ccd in metal_set])
elif ligname:
    lig_set=set([ccd for ccd in ligand_dict \
        if ligname in ligand_dict[ccd].upper()])

fasta_dict=dict()
if outfmt=='txt':
    fp=gzip.open(rootdir+"/data/protein.fasta.gz",'rt')
    for block in fp.read().split('>')[1:]:
        header,sequence=block.splitlines()
        fasta_dict[header.split()[0][1:]]=sequence
else:
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

#### parse page ####
pageLimit=200
totalNum=0
html_txt=''
fp=gzip.open(rootdir+"/data/lig_all.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    if pdbid and items[0]!=pdbid:
        continue
    if chain and items[1]!=chain and items[4]!=chain:
        continue
    if len(lig_set) and not items[3] in lig_set:
        continue
    pdb       =items[0]
    recCha    =items[1]
    bs        =items[2]
    ccd       =items[3]
    ligCha    =items[4]
    ligIdx    =items[5]
    resOrig   =items[6]
    resRenu   =items[7]
    manual    =items[8]
    moad      =items[9]
    pdbbind   =items[10]
    bindingdb =items[11]
    if baff:
        if not manual and not moad and not pdbbind and not bindingdb:
            continue
        elif baff=="manual" and not manual:
            continue
        elif baff=="moad" and not moad:
            continue
        elif baff=="pdbbind" and not pdbbind:
            continue
        elif baff=="bindingdb" and not bindingdb:
            continue            
    
    reso      =''
    csaOrig   =''
    csaRenu   =''
    ec        =''
    go        =''
    accession =''
    pmid      =''
    if not pdb+':'+recCha in chain_dict:
        if (uniprot or ecn or got or pubmed):
            continue
        totalNum+=1
    else:
        items =chain_dict[pdb+':'+recCha]
        if uniprot and items[-2]!=uniprot:
            continue
        else:
            accession=items[-2]
        if go:
            go_list=["GO:"+g for g in go.split(',')]
            go='<span title="%s">%s ...</span>'%(
                '\n'.join(go_list),go_list[0])
            if uniprot:
                go='<a href="https://ebi.ac.uk/QuickGO/annotations?geneProductId=%s" target=_blank>%s</a>'%(uniprot,go)
        else:
            go="N/A"
        if ecn and items[-4]!=ecn:
            continue
        else:
            ec=items[-4]
        if got and not go in items[-3].split(','):
            continue
        else:
            go=items[-3]
        if pubmed and items[-1]!=pubmed:
            continue
        else:
            pmid=items[-1]
        totalNum+=1
        reso   =items[0]
        if outfmt=="txt":
            csaOrig=items[1]
            csaRenu=items[2]
        else:
            if page=="last":
                if totalNum%pageLimit==1:
                    html_txt=''
            elif totalNum<=pageLimit*(int(page)-1) or pageLimit*int(page)<totalNum:
                continue
            if ec:
                ec=','.join(["<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(e,e) for e in ec.split(',')])
            else:
                ec="N/A"
            if go:
                go_list=["GO:"+g for g in go.split(',')]
                go='<span title="%s">%s ...</span>'%(
                    '\n'.join(go_list),go_list[0])
                if accession:
                    go='<a href="https://ebi.ac.uk/QuickGO/annotations?geneProductId=%s" target=_blank>%s</a>'%(accession.split(',')[0],go)
            else:
                go="N/A"
            if accession:
                accession=','.join(["<a href=https://uniprot.org/uniprot/%s target=_blank>%s</a>"%(a,a) for a in accession.split(',')])
            else:
                accession="N/A"
            if pmid:
                pmid="<a href=https://pubmed.ncbi.nlm.nih.gov/%s target=_blank>%s</a>"%(pmid,pmid)
            else:
                pmid="N/A"

    if outfmt=='txt':
        sequence  =''
        if pdb+recCha in fasta_dict:
            sequence=fasta_dict[pdb+recCha]
        html_txt+='\t'.join((pdb,recCha,reso,bs,ccd,ligCha,ligIdx,resOrig,
            resRenu,csaOrig,csaRenu,ec,go,manual,moad,pdbbind,bindingdb,
            accession,pmid,sequence))+'\n'
        continue

    name      =""
    ccd_http  =ccd
    if ccd in ligand_dict:
        name=';\n'.join(ligand_dict[ccd].split(';'))
        ccd_http='<span title="%s">%s</span>'%(name,ccd)
    
    affinity  =''
    if outfmt!='txt':
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
    bgcolor=''
    if totalNum%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    html_txt+='''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><a href="pdb.cgi?pdb=%s&chain=%s" target=_blank>%s:%s</a> (%s)</td>
    <td><span title="%s"><a href="getaid.cgi?pdb=%s&chain=%s&bs=%s" target=_blank>%s</span></td>
    <td><a href="sym.cgi?code=%s" target=_blank>%s</a></td>
    <td><a href="pdb.cgi?pdb=%s&chain=%s&idx=%s" target=_blank>%s</a></td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    totalNum,
    pdb,recCha,pdb,recCha,reso,
    resOrig,pdb,recCha,bs,bs,
    ccd,ccd_http,
    pdb,ligCha,ligIdx,ligCha,
    ec,
    go,
    accession,
    pmid,
    affinity,
    )
fp.close()

if outfmt=="txt":
    print("Content-type: text/plain\n")
    print(html_txt)
    exit()


print('''
Download all results in tab-seperated text for 
<a href="?outfmt=txt&%s" download="BioLiP.txt">%d receptor-ligand interactions</a>, whose format is explained at <a href="download/readme.txt">readme.txt</a>.<br>
Resolution -1.00 means the resolution is unavailable, e.g., for NMR structures.
Click <strong>Site #</strong> to view the binding site structure.
Hover over <strong>Site #</strong> to view the binding residues.
Hover over <strong>Ligand</strong> to view ligand details.
Hover over <strong>GO terms</strong> to view all GO terms.
<p></p>
'''%(para,totalNum))

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
<a class="hover" href="?&page=1&%s">&lt&lt</a>
<a class="hover" href="?&page=%d&%s">&lt</a>
'''%(para,page-1,para))
for p in range(page-10,page+11):
    if p<1 or p>totalPage:
        continue
    elif p==page:
        print(' %d '%(p))
    else:
        print('<a class="hover" href="?&page=%d&%s">%d</a>'%(p,para,p))
print('''
<a class="hover" href="?&page=%d&%s">&gt</a>
<a class="hover" href="?&page=last&%s">&gt&gt</a>
<form name="pform" action="qsearch.cgi">Go to page <select name="page" onchange="this.form.submit()">
'''%(page+1,para,para))
for p in range(1,totalPage+1):
    if p==page:
        print('<option value="%d" selected="selected">%d</option>'%(p,p))
    else:
        print('<option value="%d">%d</option>'%(p,p))
print('''</select>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=lig3    value='%s'>
<input type=hidden name=uniprot value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=ligname value='%s'>
<input type=hidden name=pubmed  value='%s'>
</form></center><br>'''%(pdbid,lig3,uniprot,ecn,got,ligname,pubmed))
     
print('''  
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=10% ALIGN=center><strong> PDB<br>(Resolution &#8491;) </strong></th>
    <th width=5%  ALIGN=center><strong> Site # </strong></th>
    <th width=10% ALIGN=center><strong> Ligand </strong> </th>           
    <th width=5%  ALIGN=center><strong> Ligand chain</strong> </th>           
    <th width=10% ALIGN=center><strong> EC number </strong> </th>           
    <th width=15% ALIGN=center><strong> GO terms </strong> </th>           
    <th width=10% ALIGN=center><strong> UniProt </strong> </th>           
    <th width=10% ALIGN=center><strong> PubMed </strong> </th>           
    <th width=20% ALIGN=center><strong> Binding affinity</strong> </th>           
</tr><tr ALIGN=center>
''')
print(html_txt)
print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

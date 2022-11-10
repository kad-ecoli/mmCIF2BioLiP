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
ligname=form.getfirst("ligname",'').upper().strip().strip('"').replace("'",'')
pubmed =form.getfirst("pubmed",'').strip("'")
baff   =form.getfirst("baff",'').strip().strip("'")
outfmt =form.getfirst("outfmt",'').strip().strip("'")

para_list=[]
if order:
    para_list.append("order=%s"%order)
if pdbid:
    para_list.append("pdbid=%s"%pdbid)
if chain:
    para_list.append("chain=%s"%chain)
if lig3:
    para_list.append("lig3=%s"%lig3)
elif ligname:
    para_list.append('ligname="%s"'%ligname)
if uniprot:
    para_list.append("uniprot=%s"%uniprot)
if ecn:
    para_list.append("ecn=%s"%ecn)
if got:
    para_list.append("got=%s"%got)
if baff:
    para_list.append("baff=%s"%baff)
if pubmed:
    para_list.append("pubmed=%s"%pubmed)
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
if not lig3 or not lig3 in ["rna","dna","peptide"]:
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
        if ligname in ligand_dict[ccd].replace("'",'').upper()])

fasta_dict=dict()
clust_dict=dict()
if outfmt=='txt':
    fp=gzip.open(rootdir+"/data/protein.fasta.gz",'rt')
    for block in fp.read().split('>')[1:]:
        header,sequence=block.splitlines()
        fasta_dict[header.split()[0][1:]]=sequence
else:
    if lig3 in ["rna","dna","peptide"]:
        fp=gzip.open("%s/data/%s.fasta.gz"%(rootdir,lig3),'rt')
        for block in fp.read().split('>')[1:]:
            header,sequence=block.splitlines()
            fasta_dict[header.split()[0]]=sequence
        fp.close()
        fp=gzip.open("%s/data/%s_nr.fasta.clust.gz"%(rootdir,lig3),'rt')
        for line in fp.read().splitlines():
            rep,mem=line.split('\t')
            mem_list=[rep]+mem.split(',')
            for mem in mem_list:
                clust_dict[mem]=[m.replace('_%s_'%lig3,':') for m in mem_list if m!=mem]
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
if lig3 in ["peptide","rna","dna"]:
    pageLimit=100
html_txt=''
sort_line=[]
fp=gzip.open(rootdir+"/data/lig_all.tsv.gz",'rt')
for line in fp.read().splitlines()[1:]:
#for line in stdout.decode().splitlines():
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
    manual    =items[8].replace(',',', ')
    moad      =items[9].replace(',',', ')
    pdbbind   =items[10].replace(',',', ')
    bindingdb =items[11].replace(',',', ')
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
    if got and not pdb+':'+recCha in got2chain_set:
        continue
    
    reso      =''
    csaOrig   =''
    csaRenu   =''
    ec        =''
    go        =''
    accession =''
    pmid      =''
    if not pdb+':'+recCha in chain_dict:
        if (uniprot or ecn or pubmed):
            continue
        totalNum+=1
    else:
        items =chain_dict[pdb+':'+recCha]
        if uniprot and not uniprot in items[-2].split(','):
            continue
        else:
            accession=items[-2]
        if ecn=='0' and not items[-4]:
            continue
        elif ecn and ecn!='0' and not ecn in items[-4].split(','):
            continue
        else:
            ec=items[-4]
        go=items[-3]
        if pubmed and items[-1]!=pubmed:
            continue
        else:
            pmid=items[-1]
        reso   =items[0]
        csaOrig=items[1]
        csaRenu=items[2]

    sequence  =''
    if pdb+recCha in fasta_dict:
        sequence=fasta_dict[pdb+recCha]
    items=(pdb,recCha,reso,bs,ccd,ligCha,ligIdx,resOrig,
        resRenu,csaOrig,csaRenu,ec,go,manual,moad,pdbbind,bindingdb,
        accession,pmid,sequence)
    if outfmt=='txt':
        html_txt+='\t'.join(items)+'\n'
    else:
        if order=="reso":
            sort_line.append((reso,items))
        elif order=="uniprot":
            sort_line.append((accession,items))
        elif order=="lig3":
            sort_line.append((ccd,items))
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
    
    pdb      =items[0]
    recCha   =items[1]
    reso     =items[2]
    bs       =items[3]
    ccd      =items[4]
    ligCha   =items[5]
    ligIdx   =items[6]
    resOrig  =items[7]
    resRenu  =items[8]
    csaOrig  =items[9]
    csaRenu  =items[10]
    ec       =items[11]
    go       =items[12]
    manual   =items[13]
    moad     =items[14]
    pdbbind  =items[15]
    bindingdb=items[16]
    accession=items[17]
    pmid     =items[18]
    sequence =items[19]
    
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

    name      =""
    ccd_http  =ccd
    reso="("+reso+")"
    if lig3 in ["peptide","rna","dna"]:
        ligKey='_'.join((pdb,lig3,ligCha))
        ccd_http=ccd
        if ligKey in fasta_dict:
            sequence=fasta_dict[ligKey]
            ccd_http='<br>'.join(textwrap.wrap(sequence,50))
            #ccd_http=fasta_dict[ligKey]
            if ligKey in clust_dict:
                ccd_http="&gt;"+pdb+':'+ligCha+" (identical to "+', '.join(
                    ['<a href=qsearch.cgi?lig3=%s&pdbid=%s&chain=%s target=_blank>%s</a>'%(
                    lig3,m.split(':')[0],m.split(':')[1],m
                    ) for m in clust_dict[ligKey]])+")<br>"+ccd_http
            else:
                ccd_http="&gt;"+pdb+':'+ligCha+'<br>'+ccd_http
            ccd_http='<span title="%s length=%d">%s</span>'%(
                lig3,len(sequence),ccd_http)
        reso="<br>"+reso
    else:
        if ccd in ligand_dict:
            name=';\n'.join(ligand_dict[ccd].split(';'))
            ccd_http='<span title="%s">%s</span>'%(name,ccd)
        ccd_http='<a href="sym.cgi?code=%s" target=_blank>%s</a>'%(
            ccd,ccd_http)
    
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
    title=''
    if pdb in pdb2name_dict:
        title=pdb2name_dict[pdb]
    bgcolor=''
    if l%2:
        bgcolor='BGCOLOR="#DEDEDE"'
    html_txt+='''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><span title="%s"><a href="https://rcsb.org/structure/%s" target=_blank>%s:%s</a> %s</span></td>
    <td><span title="%s"><a href="pdb.cgi?pdb=%s&chain=%s&bs=%s" target=_blank>%s</span></td>
    <td style="word-wrap: break-word">%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    l+1,
    title,pdb,pdb,recCha,reso,
    resOrig,pdb,recCha,bs,bs,
    ccd_http,
    ec,
    go,
    accession,
    pmid,
    affinity,
    )
fp.close()

print('''
Download all results in tab-seperated text for 
<a href="?outfmt=txt&%s" download="BioLiP.txt">%d receptor-ligand interactions</a>, whose format is explained at <a href="download/readme.txt">readme.txt</a>.<br>
<li>Hover over <strong>PDB</strong> to view the title of the structure.
Click <strong>PDB</strong> to view the structure at the RCSB PDB database.
Resolution -1.00 means the resolution is unavailable, e.g., for NMR structures.</li>
<li>Click <strong>Site #</strong> to view the binding site structure.
Hover over <strong>Site #</strong> to view the binding residues.</li>
'''%(para,totalNum))
if lig3 in ["peptide","rna","dna"]:
    print("<li>The <strong>sequence</strong> is converted from residues with experimentally determined coordinates in the structure; residues not observed in the 3D structure are excluded.</li>")
else:
    print("<li>Hover over <strong>Ligand</strong> to view the full ligand name. Click <strong>Ligand</strong> to view the 2D diagram and other detail information of the ligand.</li>")
print('''
<li>Hover over <strong>EC number</strong> to view the full name of enzymatic activity.</li>
<li>Hover over <strong>GO terms</strong> to view all GO terms.
<li>Hover over <strong>UniProt</strong> to view the protein name.</li>
''')


print(('''<p></p>
<form name="sform" action="qsearch.cgi">
Sort results by
<select name="order" onchange="this.form.submit()">
    <option value="pdbid">PDB ID</option>
    <option value="lig3">Ligand ID</option>
    <option value="uniprot">UniProt ID</option>
    <option value="reso">Resolution</option>
<input type=hidden name=pdbid   value='%s'>
<input type=hidden name=lig3    value='%s'>
<input type=hidden name=uniprot value='%s'>
<input type=hidden name=ecn     value='%s'>
<input type=hidden name=got     value='%s'>
<input type=hidden name=ligname value='%s'>
<input type=hidden name=pubmed  value='%s'>
<input type=hidden name=baff  value='%s'>
</form>'''%(pdbid,lig3,uniprot,ecn,got,ligname,pubmed,baff)
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


if lig3 in ["peptide","rna","dna"]:
    print('''
<style>
div.w {
  word-wrap: break-word;
}
</style>

<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=4% ALIGN=center><strong> # </strong></th>
    <th width=4% ALIGN=center><strong> PDB<br>(resolution) </strong></th>
    <th width=4% ALIGN=center><strong> Site<br># </strong></th>
    <th width=46% ALIGN=center><strong>''')
    print(lig3.upper() if lig3!="peptide" else "Peptide")
    print('''<br>sequence</strong> </th>           
    <th width=8% ALIGN=center><strong> EC<br>number </strong> </th>           
    <th width=10% ALIGN=center><strong> GO<br>terms </strong> </th>           
    <th width=8% ALIGN=center><strong> UniProt </strong> </th>           
    <th width=8% ALIGN=center><strong> PubMed </strong> </th>           
    <th width=8% ALIGN=center><strong> Binding<br>affinity</strong> </th>           
</tr><tr ALIGN=center>
''')
else:
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
print(html_txt)
print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

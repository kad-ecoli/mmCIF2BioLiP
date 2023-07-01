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

lig3=form.getfirst("code",'')
if not lig3:
    lig3=form.getfirst("lig3",'')
if lig3 in ["peptide","rna","dna","metal","regular"]:
    print('''
<style>
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
}
</style>
<table>
<tr><th>Ligand</th><th>Explanation</th></tr>
<tr><td><a href=qsearch.cgi?lig3=peptide>peptide</a></td><td>Short amino acid chain with &lt;30 residues. Old BioLiP ID: III</td></tr>
<tr><td><a href=qsearch.cgi?lig3=rna>rna</a></td><td>RNA chain. Old BioLiP ID: NUC</td></tr>
<tr><td><a href=qsearch.cgi?lig3=dna>dna</a></td><td>DNA chain. Old BioLiP ID: NUC</td></tr>
<tr><td><a href=qsearch.cgi?lig3=metal>metal</a></td><td>Metal ion</td></tr>
<tr><td><a href=qsearch.cgi?lig3=regular>regular</a></td><td>Small molecule other than a metal ion.</td></tr>
</table>
''')
elif lig3:
    txt_table='''
<style>
table, th, td {
  border: 0px solid black;
  border-collapse: collapse;
}
td {
  vertical-align: top;
  align: left;
}
</style>
<table valign=top>
'''
    freq_dict=dict()
    fp=open(rootdir+"/download/lig_frequency.txt",'r')
    for line in fp.read().splitlines()[4:]:
        items=line.split('\t')
        if len(items)==3 and items[1]==lig3:
            freq_dict[items[1]]=items[2]
    fp.close()


    filename="%s/data/smiles.tsv.gz"%rootdir
    smiles_dict=dict()
    if os.path.isfile(filename):
        fp=gzip.open(filename,'rt')
        for line in fp.read().splitlines():
            items=line.split('\t')
            if len(items)>=3 and items[0]==lig3:
                if not items[1] in smiles_dict:
                    smiles_dict[items[1]]=[items[2]]
                else:
                    smiles_dict[items[1]].append(items[2])
        fp.close()

    fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
    for line in fp.read().splitlines():
        items=line.split('\t')
        if lig3!=items[0]:
            continue
        freq="0"
        if lig3 in freq_dict and freq_dict[lig3]!='0':
            freq='<a href="qsearch.cgi?lig3=%s" target=_blank>%s</a>'%(lig3,freq_dict[lig3])
        
        SMILES=items[4]
        if len(smiles_dict):
            smiles_list=SMILES.split(';')
            SMILES="<table width=100%><tr BGCOLOR='#DEDEDE'><th>Software</th><th>SMILES</th></tr>"
            for key in smiles_list:
                key=key.strip()
                if key in smiles_dict:
                    SMILES+="<tr><td>"+'<br>'.join(smiles_dict[key])+"</td><td>"+key+"</td></tr>"
                else:
                    SMILES+="<tr><td></td><td>"+key+"</td></tr>"
            SMILES+="</table>"
        else:
            SMILES=SMILES.replace(";",";<br>")


        txt_table+='''
<tr BGCOLOR="#DEDEDE"><td width=10%><strong>PDB CCD ID: </strong></td><td width=90%><a href=https://www.rcsb.org/ligand/$ccd target=_blank>$ccd</a></td></tr>
<tr><td><strong>Number of entries in BioLiP: </strong></td><td>$freq</td></tr>
<tr BGCOLOR="#DEDEDE"><td><strong>Chemical formula: </strong></td><td>$formula</td></tr>
<tr><td><strong>InChI: </strong></td><td>$InChI</td></tr>
<tr BGCOLOR="#DEDEDE"><td><strong>InChIKey: </strong></td><td>$InChIKey</td></tr>
<tr><td><strong>SMILES: </strong></td><td>$SMILES</td></tr>
<tr BGCOLOR="#DEDEDE"><td><strong>Name:</strong></td><td>$name</td></tr>
      '''.replace("$ccd",lig3
        ).replace("$freq",freq
        ).replace("$formula",items[1]
        ).replace("$InChIKey",items[3]
        ).replace("$InChI",items[2]
        ).replace("$SMILES",SMILES
        ).replace("$name",items[5].replace(';',';<br>'))
        BGCOLOR=""
        if items[6]:
            txt_table+='''<tr><td><strong>ChEMBL: </strong></td><td><a href="https://www.ebi.ac.uk/chembl/compound_report_card/$ChEMBL" target=_blank>$ChEMBL</a></td></tr>
            '''.replace("$ChEMBL",items[6])
            BGCOLOR='BGCOLOR="#DEDEDE"'
        if items[7]:
            txt_table+='''<tr $BGCOLOR><td><strong>DrugBank: </strong></td><td><a href="https://go.drugbank.com/drugs/$DrugBank" target=_blank>$DrugBank</a></td></tr>
            '''.replace("$DrugBank",items[7]
              ).replace("$BGCOLOR", BGCOLOR)
            if BGCOLOR:
                BGCOLOR=""
            else:
                BGCOLOR='BGCOLOR="#DEDEDE"'
        if items[8]:
            txt_table+='''<tr $BGCOLOR><td><strong>ZINC: </strong></td><td><a href="https://zinc.docking.org/substances/$ZINC" target=_blank>$ZINC</a></td></tr>
            '''.replace("$ZINC",items[8]
              ).replace("$BGCOLOR", BGCOLOR)
        svg="https://cdn.rcsb.org/images/ccd/labeled/%s/%s.svg"%(lig3[0],lig3)
        print("<p><ul><a href=%s target=_blank><img src=%s alt='' width=300></a><br>"%(svg,svg))
        print("View %s at the <a href=https://rcsb.org/ligand/%s target=_blank>PDB</a> and <a href=qsearch.cgi?lig3=%s>BioLiP</a> database</ul></p>"%(
            lig3,lig3,lig3))
    fp.close()
    print(txt_table+"</table>")
else:
    fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
    ligand_list=[]
    for line in fp.read().splitlines()[1:]:
        ligand_list.append(line.split('\t')[0])
    fp.close()
    import random
    lig3=random.choice(ligand_list)
    
    print('''No ligand code provided. You may <a href="%s?code=%s">[browse a random ligand]</a>
<meta http-equiv="refresh" content="3; url='%s?code=%s'" /><br>
This page will be redicted in 3 seconds.
'''%(os.path.basename(__file__),lig3,
     os.path.basename(__file__),lig3))

if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

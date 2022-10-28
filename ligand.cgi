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

form    = cgi.FieldStorage()
lig3    = form.getfirst("code",'').upper()
if not lig3:
    lig3= form.getfirst("lig3",'').upper()
formula = form.getfirst("formula",'').upper()
inchi   = form.getfirst("inchi",'').upper()
if inchi and not inchi.startswith("INCHI="):
    inchi="INCHI="+inchi
inchikey= form.getfirst("inchikey",'').upper()
smiles  = form.getfirst("smiles",'').upper()
ligname = form.getfirst("ligname",'').upper()

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
<h4>Search Ligand</h4>
<FORM id="form1" name="form1" METHOD="POST" ENCTYPE="MULTIPART/FORM-DATA" ACTION="ligand.cgi">
    <span title="Example: 5GP">
    Ligand ID
    <input type="text" name="code" value="" placeholder="%s" size=3>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: C10 H14 N5 O8 P">
    Formula
    <input type="text" name="formula" value="" placeholder="%s" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: InChI=1S/C10H14N5O8P/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(23-9)1-22-24(19,20)21/h2-3,5-6,9,16-17H,1H2,(H2,19,20,21)(H3,11,13,14,18)/t3-,5-,6-,9-/m1/s1">
    InChI
    <input type="text" name="inchi" value="" placeholder="%s" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: RQFCJASXJCIDSX-UUOKFMHZSA-N">
    InChIKey
    <input type="text" name="inchikey" value="" placeholder="%s" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N">
    SMILES
    <input type="text" name="smiles" value="" placeholder="%s" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: GUANOSINE-5'-MONOPHOSPHATE">
    Ligand name
    <input type="text" name="ligname" value="" placeholder="%s" size=20>
    </span>
    &nbsp;&nbsp;
    
    <INPUT TYPE="submit" VALUE="Submit">
    &nbsp;
    <INPUT TYPE="reset" VALUE="Clear">
</FORM>
<p></p>
'''%(lig3,formula,inchi,inchikey,smiles,ligname))

print('''
<h4>Search Ligand Result</h4>
Click the corresponding <strong> Ligand ID</strong> to search BioLiP. The ligand ID follows the <a href="https://www.wwpdb.org/data/ccd" target=_blank>Chemical Component Dictionary (CCD)</a> used by the PDB database.<br>
Click the corresponding <strong> Ligand Name</strong> to visualize the ligand and display its names/synonyms.<br>
If multiple SMILES strings exists for the same ligand, different SMILES are separated by semicolon ";".
<p></p>
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=4% ALIGN=center><strong> # </srong></th>
    <th width=6% ALIGN=center><strong> Ligand ID </srong></th>
    <th width=10% ALIGN=center><strong> Formula </srong></th>
    <th width=10% ALIGN=center><strong> InChI </srong></th>
    <th width=10% ALIGN=center><strong> InChIKey </srong></th>
    <th width=10% ALIGN=center><strong> SMILES </srong></th>
    <th width=50% ALIGN=left>  <strong> Ligand Name </srong> </th>           
</tr><tr ALIGN=center>
''')

count=0
fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
formula_query_set=set(formula.split())
for line in fp.read().splitlines()[1:]:
    items=line.split('\t')
    if lig3 and items[0]!=lig3:
        continue
    if formula:
        formula_template_set=set(items[1].split())
        if len(formula_query_set)!=len(formula_template_set):
            continue
        if len(formula_query_set.intersection(formula_template_set)
            )!=len(formula_template_set):
            continue
    if inchi and items[2].upper()!=inchi:
        continue
    if inchikey and items[3]!=inchikey:
        continue
    if smiles:
        smiles_list=items[4].split(';')
        if sum([s.strip().upper()==smiles for s in smiles_list])==0:
            continue
    if ligname and not ligname in items[5].upper():
        continue
    count+=1
    
    inchi_list=[]
    for i in range(0,len(items[2]),25):
        start=i
        end  =i+25
        inchi_list.append(items[2][start:end])
    items[2]='<br>'.join(inchi_list)
    
    smiles_list=[]
    for i in range(0,len(items[4]),40):
        start=i
        end  =i+40
        smiles_list.append(items[4][start:end])
    items[4]='<br>'.join(smiles_list)
    items[4].replace('; ',';<br>')
    
    bgcolor=''
    if count%2==0:
        bgcolor='BGCOLOR="#DEDEDE"'
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td><a href="qsearch.cgi?lig3=%s" target="_blank">%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td>%s</td>
    <td ALIGN=left><a href="sym.cgi?code=%s" target="_blank">%s</td>
</tr>
'''%(bgcolor,
    count,
    items[0],items[0],
    items[1],
    items[2],
    items[3],
    items[4],
    items[0],items[5].replace(';',';<br>')
    ))
fp.close()


print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

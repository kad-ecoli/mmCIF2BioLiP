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
freq_dict=dict()
fp=open(rootdir+"/download/lig_frequency.txt")
for line in fp.read().splitlines()[4:]:
    items=line.split('\t')
    ccd=items[1]
    freq=items[2]
    freq_dict[ccd]=freq

fp=gzip.open(rootdir+"/data/ligand.tsv.gz",'rt')
lines=fp.read().splitlines()[1:]
fp.close()

totalNum=len(lines)

print('''
<h4>Search Ligand</h4>
<FORM id="form1" name="form1" METHOD="POST" ENCTYPE="MULTIPART/FORM-DATA" ACTION="ligand.cgi">
    <span title="Example: 5GP">
    Ligand ID
    <input type="text" name="code" value="" size=3>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: C10 H14 N5 O8 P">
    Formula
    <input type="text" name="formula" value="" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: InChI=1S/C10H14N5O8P/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(23-9)1-22-24(19,20)21/h2-3,5-6,9,16-17H,1H2,(H2,19,20,21)(H3,11,13,14,18)/t3-,5-,6-,9-/m1/s1">
    InChI
    <input type="text" name="inchi" value="" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: RQFCJASXJCIDSX-UUOKFMHZSA-N">
    InChIKey
    <input type="text" name="inchikey" value="" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N">
    SMILES
    <input type="text" name="smiles" value="" size=10>
    </span>
    &nbsp;&nbsp;
    
    <span title="Example: GUANOSINE-5'-MONOPHOSPHATE">
    Ligand name
    <input type="text" name="ligname" value="" size=20>
    </span>
    &nbsp;&nbsp;
    
    <INPUT TYPE="submit" VALUE="Submit">
    &nbsp;
    <INPUT TYPE="reset" VALUE="Clear">
</FORM>
<p></p>

<h4>Browse Ligand</h4>
<strong> %d </strong> ligands in BioLiP.<br>
Click the corresponding <strong> Ligand ID</strong> to search BioLiP. The ligand ID follows the <a href="https://www.wwpdb.org/data/ccd" target=_blank>Chemical Component Dictionary (CCD)</a> used by the PDB database.<br>
<strong>Count</strong> refers to the number of BioLiP entries associated with the ligand. The full statistics is available at <a href=download/lig_frequency.txt>lig_frequency.txt</a><br>
Click the corresponding <strong> Ligand Name</strong> to visualize the ligand and display its names/synonyms.<br>
<p></p>
'''%totalNum)

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
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=8% ALIGN=center><strong> Ligand ID </strong></th>
    <th width=7% ALIGN=center><strong> Count </strong></th>
    <th width=80% ALIGN=left>  <strong> Ligand Name </strong> </th>           
</tr><tr ALIGN=center>
''')
for l in range(pageLimit*(page-1),pageLimit*page+1):
    if l>=totalNum:
        continue
    items=lines[l].split('\t')
    ccd  =items[0]
    freq ='0'
    if ccd in freq_dict:
        freq = freq_dict[ccd]
    name =items[-1]
    name =';<br>'.join(name.split(';'))
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
'''%(bgcolor,l+1,ccd,ccd,freq,ccd,name))


print("</table>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import subprocess
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

def ExitWithError(msg,html_footer):
    print("ERROR!")
    print(msg)
    print("<p></p><a href=.>[Back]</a>")
    if len(html_footer):
        print(html_footer)
    else:
        print("</body> </html>")
    exit()

#### read cgi parameters ####

form = cgi.FieldStorage()
sequence=form.getfirst("sequence",'').strip()
if not sequence:
    sequence=form.getfirst("seq_file",'').strip()
seq_type =form.getfirst("seq_type",'').strip().lower()

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
if not seq_type in ["protein","rna","dna","peptide"]:
    ExitWithError("Sequence type must be one of the following: protein, rna, dna, peptide",html_footer)
header=seq_type
txt=''
for line in sequence.splitlines():
    line=line.strip()
    if line.startswith('>'):
        if header[0]=='>':
            print("ERROR! only one sequence allowed per search")
            exit()
        else:
            header=line
    else:
        txt+=line.upper()
if header[0]=='>':
    header=header[1:]
if seq_type in ["rna","dna"]:
    sequence=txt.lower()
else:
    sequence=txt.upper()
if len(set(txt).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))):
    ExitWithError("Unknown residue type "+' '.join(set(txt
        ).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))),html_footer)
if len(sequence)>1500:
    ExitWithError("Unable to handle sequence with %d &gt; 1500 residues"%len(sequence),html_footer)
seqID=100
if seq_type=="protein":
    seqID=90
print("Search sequence (length=%d) through non-redundant sequence database <a href=data/%s_nr.fasta.gz>%s_nr.fasta.gz</a> clustered at %d%% identity cutoff.<br>"%(len(sequence),seq_type,seq_type,seqID))
print('&gt;'+header+'<br>')
print('<br>'.join(textwrap.wrap(sequence,80))+'<p></p>')

blast="blastp"
if seq_type.endswith("na"):
    blast="blastn"
cmd="echo %s|%s/script/%s -db %s/data/%s_nr -max_target_seqs 1000 -outfmt '6 sacc slen evalue nident length' "%(sequence,rootdir,blast,rootdir,seq_type)
score_name="E-value"
if len(sequence)<30:
    cmd="echo %s | %s/script/NWalign - %s/data/%s_nr.fasta.gz|sort -k3nr|head -1000"%(sequence,rootdir,rootdir,seq_type)
    score_name="Alignment<br>score"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
lines=stdout.decode().splitlines()
print('''
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=10% ALIGN=center><strong> Hit </strong></th>
    <th width=10% ALIGN=center><strong> Hit<br>length </strong></th>
    <th width=10% ALIGN=center><strong> Aligned<br>length </strong></th>
    <th width=10% ALIGN=center><strong> Identity<br>(normalized by query)</strong> </th>           
    <th width=10% ALIGN=center><strong> Identity<br>(normalized by hit)</strong> </th>           
    <th width=10% ALIGN=center><strong> Identity (normalized<br>by aligned length)</strong> </th>           
    <th width=10% ALIGN=center><strong> '''+score_name+ '''</strong> </th>           
    <th width=25% ALIGN=center><strong> Homologs<br>to hit</strong> </th>           
</tr><tr ALIGN=center>
''')
hit2chain_dict=dict()
hit2clust_dict=dict()
if len(lines):
    if seq_type=="protein":
        cmd="zcat %s/data/%s.fasta.gz|grep '>'|cut -f1,2|sed 's/^>//g'"%(rootdir,seq_type)
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        for line in stdout.decode().splitlines():
            hit,chainID=line.split('\t')
            pdbid=hit[:-len(chainID)]
            hit2chain_dict[hit]=(pdbid,chainID)
    fp=gzip.open("%s/data/%s_nr.fasta.clust.gz"%(rootdir,seq_type),'rt')
    for line in fp.read().splitlines():
        rep,mem=line.split('\t')
        hit2clust_dict[rep]=mem.split(',')
    fp.close()

totalNum=0
sacc_list=[]
for line in lines:
    sacc,slen,evalue,nident,Lali=line.split('\t')
    if sacc in sacc_list:
        continue
    totalNum+=1
    bgcolor=''
    if totalNum%2==0:
        bgcolor='BGCOLOR="#DEDEDE"'
    slen=int(slen)
    nident=float(nident)
    Lali=int(Lali)

    
    if seq_type=="protein":
        pdbid,chainID=hit2chain_dict[sacc]
        hit="<a href=qsearch.cgi?pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
            pdbid,chainID,pdbid,chainID)
    else:
        pdbid,chainID=sacc.split('_%s_'%seq_type)
        hit="<a href=qsearch.cgi?lig3=%s&pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
            seq_type,pdbid,chainID,pdbid,chainID)
    homolog_list=[]
    if sacc in hit2clust_dict:
        for mem in hit2clust_dict[sacc]:
            if seq_type=="protein":
                pdbid,chainID=hit2chain_dict[mem]
                homolog_list.append("<a href=qsearch.cgi?pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
                    pdbid,chainID,pdbid,chainID))
            else:
                pdbid,chainID=mem.split('_%s_'%seq_type)
                homolog_list.append("<a href=qsearch.cgi?lig3=%s&pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
                    seq_type,pdbid,chainID,pdbid,chainID))
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td>%s</td>
    <td>%d</td>
    <td>%d</td>
    <td>%.4f</td>
    <td>%.4f</td>
    <td>%.4f</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    totalNum,
    hit,
    slen,
    Lali,
    nident/len(sequence),
    nident/slen, 
    nident/Lali,
    evalue, 
    ', '.join(homolog_list)))
print("</table><p></p><a href=.>[Back]</a>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

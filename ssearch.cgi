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
    print("<p></p><a href=ssearch.html>[Back to search]</a>")
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
sequence=txt
if len(set(txt).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))):
    ExitWithError("Unknown residue type "+' '.join(set(txt
        ).difference(set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))),html_footer)
if len(sequence)>1500:
    ExitWithError("Unable to handle sequence with %d &gt; 1500 residues"%len(sequence),html_footer)
print("Search sequence through non-redundant sequence database <a href=data/%s_nr.fasta.gz>%s_nr.fasta.gz</a><br>"%(seq_type,seq_type))
print('&gt;'+header+'<br>')
print('<br>'.join(textwrap.wrap(sequence,80))+'<p></p>')

blast="blastp"
if seq_type.endswith("na"):
    blast="blastn"
cmd="echo %s|%s/script/%s -db %s/data/%s_nr -outfmt '6 sacc slen evalue nident' "%(sequence,rootdir,blast,rootdir,seq_type)
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
stdout=stdout.decode()
print('<br>'.join(stdout.splitlines()))
    

print("<p></p><a href=ssearch.html>[Back to search]</a>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")

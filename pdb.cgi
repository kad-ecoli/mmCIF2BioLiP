#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import subprocess
import textwrap
import tarfile
import gzip

rootdir=os.path.dirname(os.path.abspath(__file__))

def display_ligand(pdbid,asym_id,lig3,ligIdx):
    reso  =''
    pubmed=''
    cmd="zcat %s/data/pdb_all.tsv.gz |cut -f1,3,9|uniq|grep -P '^%s\\t'"%(
        rootdir,pdbid)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    items=stdout.decode().split('\t')
    if len(items)>=3:
        reso,pubmed=items[1:3]
        if reso=="-1.00":
            reso="N/A"
        else:
            reso=reso+" &#8491;"
    
    prefix='_'.join((pdbid,lig3,asym_id,ligIdx))
    filename="%s/output/%s.pdb.gz"%(rootdir,prefix)
    if not os.path.isfile(filename):
        divided=pdbid[-3:-1]
        tar = tarfile.open("%s/weekly/ligand_%s.tar.bz2"%(rootdir,divided))
        fin=tar.extractfile("ligand/%s.pdb"%prefix)
        lines=fin.read().decode().splitlines()
        txt=''
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" L"+line[22:]+'\n'
            else:
                txt+=line+'\n'
        fout=gzip.open(filename,'wt')
        fout.write(txt)
        fout.close()
        fin.close()
        tar.close()
    script=''
    if lig3 in ["peptide"]:
        script="cartoons; color group;"
    elif lig3 in ["rna","dna"]:
        script="cartoons; color group; spacefill off; wireframe off;"

    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">3D structure</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>PDB</strong></td><td><a href=qsearch.cgi?pdbid=$pdbid target=_blank>$pdbid</a></td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><strong>Chain</strong></td><td><a href=qsearch.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$asym_id</a></td></tr>
    <tr><td align=center><strong>Resolution</strong></td><td>$reso</a></td></tr>
    <tr BGCOLOR="#DEDEDE" align=center><td align=center><strong>3D<br>structure</strong></td><td>
    <table><tr><td>

<script type="text/javascript"> 
$(document).ready(function()
{
    Info = {
        width: 400,
        height: 400,
        j2sPath: "jsmol/j2s",
        script: "load output/$prefix.pdb.gz; color background black; $script"
    }
    $("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))
});
</script>
<span id=mydiv></span>
    </td><td align=left>
[<a href="javascript:Jmol.script(jmolApplet0, 'spin on')">Spin on</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'spin off')">Spin off</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'Reset')">Reset orientation</a>]<p></p>
[<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay true')">High quality</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay false')">Low quality</a>]<p></p>
[<a href="javascript:Jmol.script(jmolApplet0, 'color background white')">White quality</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'color background black')">Black quality</a>]<p></p>
[<a href=output/$prefix.pdb.gz>Download</a>]
    </td></tr></table>


    </td></tr>
    </table>
</div>
</td></tr>
'''.replace("$pdbid",pdbid
  ).replace("$asym_id",asym_id
  ).replace("$reso",reso
  ).replace("$prefix",prefix
  ).replace("$script",script
  ))


    cmd="zcat %s/data/lig_all.tsv.gz|grep -P '^%s\\t\w+\\tBS\d+\\t%s\\t%s\\t%s\\t'"%(
        rootdir,pdbid,lig3,asym_id,ligIdx)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    lines=stdout.decode().splitlines()
    if len(lines):
        print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Interaction with protein receptor</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr BGCOLOR="#DEDEDE" align=center>
        <th width=10%><strong>Receptor chain</strong></th>
        <th width=30%><strong>Binding residues on receptor<br>(original residue number in PDB)</strong></th>
        <th width=30%><strong>Binding residues on receptor<br>(residue number reindexed from 1)</strong></th>
        <th width=30%><strong>Binding affinity</strong></th>
    </tr>
''')
        for l,line in enumerate(lines):
            items    =line.split('\t')
            recCha   =items[1]
            bs       =items[2]
            resOrig  =items[6]
            resRenu  =items[7]
            manual   =items[8]
            moad     =items[9]
            pdbbind  =items[10]
            bindingdb=items[11]

            baff_list=[]
            if manual:
                baff_list.append("Manual survey: "+manual)
            if moad:
                baff_list.append("<a href=http://bindingmoad.org/pdbrecords/index/%s target=_blank>MOAD</a>: %s"%(pdbid,moad))
            if pdbbind:
                baff_list.append("<a href=http://pdbbind.org.cn/quickpdb.php?quickpdb=%s target=_blank>PDBbind-CN</a>: %s"%(pdbid,pdbbind))
            if bindingdb:
                baff_list.append("BindingDB: "+bindingdb)

            bgcolor=''
            if l%2==1:
                bgcolor=' BGCOLOR="#DEDEDE" '
            print('''
    <tr $bgcolor align=center>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$pdbid:$recCha</a></span></td>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$resOrig</a></span></td>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$resRenu</a></span></td>
        <td>$baff</td>
    </tr>
            '''.replace("$bgcolor",bgcolor
              ).replace("$pdbid",pdbid
              ).replace("$recCha",recCha
              ).replace("$bs",bs
              ).replace("$resOrig",resOrig
              ).replace("$resRenu",resRenu
              ).replace("$baff",'<br>'.join(baff_list)))

        print('''   </table>
</div>
</td></tr>
''')

    return pubmed

    

def display_polymer_ligand(pdbid,asym_id,lig3):
    code=lig3.upper()
    if lig3=="peptide":
        code="peptide"
    prefix="%s_%s_%s"%(pdbid,lig3,asym_id)
    cmd="zcat %s/data/%s.fasta.gz|grep -PA1 '^>%s\\t'|tail -1"%(
        rootdir,lig3,prefix)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    sequence=stdout.decode().strip()
    seq_txt='<br>'.join(textwrap.wrap(sequence,50))

    cmd="zcat %s/data/%s_nr.fasta.clust.gz|sed 's/\\t/,/g'|grep -P '\\b%s\\b'"%(
        rootdir,lig3,prefix)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    homolog_link=''
    for homolog in stdout.decode().strip().split(','):
        if homolog==prefix or not homolog.strip():
            continue
        homo_pdbid,homo_asym_id=homolog.split('_'+lig3+'_')
        homolog_link+=", <a href=pdb.cgi?pdb=%s&chain=%s&lig3=%s&idx=0 target=_blank>%s:%s</a>"%(
            homo_pdbid,homo_asym_id,lig3,homo_pdbid,homo_asym_id)
    if homolog_link:
        homolog_link='<tr BGCOLOR="#DEDEDE"><td>(Identical to '+ \
        homolog_link[1:]+")</td></tr>"

    print('''
<tr><td><h1 align=center>Structure of PDB $pdbid Chain $asym_id</h1></td></tr>
<tr><td>
<div id="headerDiv">
    <div id="titleText">Sequence</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td>&gt;$prefix (length=$L) [<a href=ssearch.cgi?seq_type=$lig3&sequence=$sequence>Search $code sequence</a>]</td></tr>
    <tr><td><span title="Only residues with experimentally determined coordinates are included. Residues unobserved in the structure are excluded.">$seq_txt</span></td></tr>
    $homolog_link
    </table>
</div>
</td></tr>
'''.replace("$pdbid",pdbid
  ).replace("$asym_id",asym_id
  ).replace("$prefix",prefix
  ).replace("$L",str(len(sequence))
  ).replace("$lig3",lig3
  ).replace("$code",code
  ).replace("$sequence",sequence
  ).replace("$seq_txt",seq_txt
  ).replace("$homolog_link",homolog_link
  ))
    return display_ligand(pdbid,asym_id,lig3,'0')

def display_regular_ligand(pdbid,asym_id,lig3,ligIdx):
    if not ligIdx:
        ligIdx='1'
    lig3=lig3.upper()

    formula=''
    InChI=''
    InChIKey=''
    SMILES=''
    name=''
    cmd="zcat %s/data/ligand.tsv.gz|grep -P '^%s\\t'"%(rootdir,lig3)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    items=stdout.decode().split('\t')
    if len(items)>=6:
        formula,InChI,InChIKey,SMILES,name=items[1:6]
    svg="https://cdn.rcsb.org/images/ccd/labeled/%s/%s.svg"%(lig3[0],lig3)
    print('''
<tr><td><h1 align=center>Structure of PDB $pdbid Chain $asym_id ligand $lig3</h1></td></tr>
<tr><td>
<div id="headerDiv">
    <div id="titleText">Chemical information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr align=center><td width=10%><strong>2D<br>diagram</strong></td><td><a href=$svg target=_blank><img src=$svg width=400></a></td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><a href=https://wwpdb.org/data/ccd target=_blank><strong>Ligand ID</strong></a></td><td><span title="The ligand ID follows the
Chemical Component Dictionary (CCD)
used by the PDB database."><a href=https://rcsb.org/ligand/$lig3>$lig3</a></span></td></tr>
    <tr><td align=center><a href=https://inchi-trust.org target=_blank><strong>InChI</strong></a></td><td>$InChI</td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><a href=https://inchi-trust.org target=_blank><strong>InChIKey</strong></a></td><td>$InChIKey</td></tr>
    <tr><td align=center><strong>SMILES</strong></td><td>$SMILES</td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center><strong>Formula</strong></td><td>$formula</td></tr>
    <tr><td align=center><strong>Name</strong></td><td>$name</td></tr>
    </table>
</div>
</td></tr>
'''.replace("$pdbid",pdbid
  ).replace("$asym_id",asym_id
  ).replace("$lig3",lig3
  ).replace("$svg",svg
  ).replace("$InChIKey",InChIKey
  ).replace("$InChI",InChI
  ).replace("$SMILES",SMILES.replace(';',';<br>')
  ).replace("$formula",formula
  ).replace("$name",name.replace(';',';<br>')
  ))
    return display_ligand(pdbid,asym_id,lig3,ligIdx)

def display_ec(ec,csaOrig,csaRenu):
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Enzymatic activity</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
''')
    if csaOrig:
        print('''
    <tr align=center>
        <td width=10%><strong>Catalytic site (original residue number in PDB)</strong></td>
        <td width=90%><a href="https://www.ebi.ac.uk/thornton-srv/m-csa/search/?s=$pdbid" target=_blank>$csaOrig</a></td>
    </tr>
    <tr BGCOLOR="#DEDEDE" align=center>
        <td width=10%><strong>Catalytic site (residue number reindexed from 1)</strong></td>
        <td width=90%><a href="https://www.ebi.ac.uk/thornton-srv/m-csa/search/?s=$pdbid" target=_blank>$csaRenu</a></td>
    </tr>
        '''.replace("$csaOrig",csaOrig
          ).replace("$csaRenu",csaRenu
          ).replace("$pdbid",pdbid))
    if ec:
        ec_list=["<a href=https://enzyme.expasy.org/EC/%s target=_blank>%s</a>"%(
            e,e) for e in ec.replace(' ','').split(',')]
        print('''
    <tr align=center>
        <td width=10%><strong>Enzyme Commision number</strong></td>
        <td width=90%>$ec</td>
    </tr>
        '''.replace("$ec",'; '.join(ec_list)))
    print('''   </table>
</div>
</td></tr>
''')

def display_go(go):
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Gene Ontology</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr BGCOLOR="#DEDEDE" align=center>
        <th width=10%><strong>GO term</strong></th>
        <th width=80%><strong>Name</strong></th>
        <th width=10%><strong>Aspect<br>chain</strong></th>
    </tr>
''')
    for l,term in enumerate(go.split(',')):
        term="GO:"+term.strip()
        bgcolor=''
        if l%2==1:
            bgcolor=' BGCOLOR="#DEDEDE" '
        print('''
    <tr $bgcolor align=center>
        <td><a href="https://ebi.ac.uk/QuickGO/term/$term" target=_blank>$term</td>
        <td>$name</td>
        <td>$aspect</td>
    </tr>
        '''.replace("$bgcolor",bgcolor
          ).replace("$term",term))

    print('''   </table>
</div>
</td></tr>
''')

def display_protein_receptor(pdbid,asym_id):
    prefix="%s%s"%(pdbid,asym_id)
    cmd="zcat %s/data/protein.fasta.gz|grep -PA1 '^>%s\\t'|tail -1"%(
        rootdir,prefix)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    sequence=stdout.decode().strip()
    seq_txt='<br>'.join(textwrap.wrap(sequence,50))
    print('''
<tr><td><h1 align=center>Structure of PDB $pdbid Chain $asym_id</h1></td></tr>
<tr><td>
<div id="headerDiv">
    <div id="titleText">Receptor sequence</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td>&gt;$prefix (length=$L) [<a href=ssearch.cgi?seq_type=protein&sequence=$sequence target=_blank>Search protein sequence</a>]</td></tr>
    <tr><td><span title="Only residues with experimentally determined coordinates are included. Residues unobserved in the structure are excluded.">$seq_txt</span></td></tr>
    </table>
</div>
</td></tr>
'''.replace("$pdbid",pdbid
  ).replace("$asym_id",asym_id
  ).replace("$prefix",prefix
  ).replace("$L",str(len(sequence))
  ).replace("$sequence",sequence
  ).replace("$seq_txt",seq_txt
  ))

    reso   =''
    csaOrig=''
    csaRenu=''
    ec     =''
    go     =''
    uniprot=''
    pubmed =''
    cmd="zcat %s/data/pdb_all.tsv.gz |grep -P '^%s\\t%s\\t'|head -1"%(
        rootdir,pdbid,asym_id)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    items=stdout.decode().split('\t')
    if len(items)>=9:
        reso,csaOrig,csaRenu,ec,go,uniprot,pubmed=items[2:9]
        if reso=="-1.00":
            reso="N/A"
        else:
            reso=reso+" &#8491;"
    filename="%s/output/%s.pdb.gz"%(rootdir,prefix)
    if not os.path.isfile(filename):
        divided=pdbid[-3:-1]
        tar = tarfile.open("%s/weekly/receptor_%s.tar.bz2"%(rootdir,divided))
        fin=tar.extractfile("receptor/%s.pdb"%prefix)
        lines=fin.read().decode().splitlines()
        txt=''
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" R"+line[22:]+'\n'
            else:
                txt+=line+'\n'
        fout=gzip.open(filename,'wt')
        fout.write(txt)
        fout.close()
        fin.close()
        tar.close()

    script=''
    explainLabel=''
    if csaOrig:
        explainLabel='Catalytic site residues are labeled in the structure<br>'
        resi_list=[r[1:] for r in csaOrig.split()]
        script="select "+','.join(resi_list)+"; spacefill 25%; wireframe 50; color group;"
        for resi in resi_list:
            script+="select "+resi+" and *.ca; label %m%R; color label magenta;"
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">3D structure</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%>PDB</td><td><a href=qsearch.cgi?pdbid=$pdbid target=_blank>$pdbid</a></td></tr>
    <tr BGCOLOR="#DEDEDE"><td align=center>Chain</td><td><a href=qsearch.cgi?pdbid=$pdbid&chain=$asym_id target=_blank>$asym_id</a></td></tr>
    <tr><td align=center>Resolution</td><td>$reso</a></td></tr>
    <tr BGCOLOR="#DEDEDE" align=center><td align=center>3D<br>structure</td><td>
    <table><tr><td>

<script type="text/javascript"> 
$(document).ready(function()
{
    Info = {
        width: 400,
        height: 400,
        j2sPath: "jsmol/j2s",
        script: "load output/$prefix.pdb.gz; color background black; cartoons; color group; spacefill off; wireframe off; $script;"
    }
    $("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))
});
</script>
<span id=mydiv></span>
    </td><td align=left>
$explainLabel
[<a href="javascript:Jmol.script(jmolApplet0, 'spin on')">Spin on</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'spin off')">Spin off</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'Reset')">Reset orientation</a>]<p></p>
[<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay true')">High quality</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay false')">Low quality</a>]<p></p>
[<a href="javascript:Jmol.script(jmolApplet0, 'color background white')">White quality</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'color background black')">Black quality</a>]<p></p>
[<a href=output/$prefix.pdb.gz>Download</a>]
    </td></tr></table>


    </td></tr>
    </table>
</div>
</td></tr>
'''.replace("$pdbid",pdbid
  ).replace("$asym_id",asym_id
  ).replace("$reso",reso
  ).replace("$prefix",prefix
  ).replace("$explainLabel",explainLabel
  ).replace("$script",script
  ))

    if ec:
        display_ec(ec,csaOrig,csaRenu)

    cmd="zcat %s/data/lig_all.tsv.gz|grep -P '^%s\\t%s\\t'"%(
        rootdir,pdbid,asym_id)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    lines=stdout.decode().splitlines()
    if len(lines):
        print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Interaction with ligand</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr BGCOLOR="#DEDEDE" align=center>
        <th width=5%><strong>Site<br>#</strong></th>
        <th width=5%><strong>Ligand</strong></th>
        <th width=5%><strong>Ligand<br>chain</strong></th>
        <th width=29%><strong>Binding residues on receptor<br>(original residue number in PDB)</strong></th>
        <th width=29%><strong>Binding residues on receptor<br>(residue number reindexed from 1)</strong></th>
        <th width=27%><strong>Binding affinity</strong></th>
    </tr>
''')
        for l,line in enumerate(lines):
            items    =line.split('\t')
            recCha   =items[1]
            bs       =items[2]
            ccd      =items[3]
            ligCha   =items[4]
            ligIdx   =items[5]
            resOrig  =items[6]
            resRenu  =items[7]
            manual   =items[8]
            moad     =items[9]
            pdbbind  =items[10]
            bindingdb=items[11]

            baff_list=[]
            if manual:
                baff_list.append("Manual survey: "+manual)
            if moad:
                baff_list.append("<a href=http://bindingmoad.org/pdbrecords/index/%s target=_blank>MOAD</a>: %s"%(pdbid,moad))
            if pdbbind:
                baff_list.append("<a href=http://pdbbind.org.cn/quickpdb.php?quickpdb=%s target=_blank>PDBbind-CN</a>: %s"%(pdbid,pdbbind))
            if bindingdb:
                baff_list.append("BindingDB: "+bindingdb)

            bgcolor=''
            if l%2==1:
                bgcolor=' BGCOLOR="#DEDEDE" '
            print('''
    <tr $bgcolor align=center>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$bs</a></span></td>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$ccd</a></span></td>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$ligCha</a></span></td>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$resOrig</a></span></td>
        <td><span title="Click to view binding site"><a href="pdb.cgi?pdb=$pdbid&chain=$recCha&bs=$bs" target=_blank>$resRenu</a></span></td>
        <td>$baff</td>
    </tr>
            '''.replace("$bgcolor",bgcolor
              ).replace("$pdbid",pdbid
              ).replace("$recCha",recCha
              ).replace("$ccd",ccd
              ).replace("$ligCha",recCha
              ).replace("$bs",bs
              ).replace("$resOrig",resOrig
              ).replace("$resRenu",resRenu
              ).replace("$baff",'<br>'.join(baff_list)))

        print('''   </table>
</div>
</td></tr>
''')

    if go:
        display_go(go)
    return pubmed,uniprot

def display_interaction(pdbid,asym_id,bs):    
    prefix="%s%s"%(pdbid,asym_id)
    cmd="zcat %s/data/protein.fasta.gz|grep -PA1 '^>%s\\t'|tail -1"%(
        rootdir,prefix)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    sequence=stdout.decode().strip()
    seq_txt='<br>'.join(textwrap.wrap(sequence,50))
    print('''
<tr><td><h1 align=center>Structure of PDB $pdbid Chain $asym_id Binding Site $bs</h1></td></tr>
<tr><td>
<div id="headerDiv">
    <div id="titleText">Receptor Information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td>&gt;$pdbid Chain $asym_id (length=$L)
    [<a href=ssearch.cgi?seq_type=protein&sequence=$sequence target=_blank>Search protein sequence</a>]
    [<a href=output/$prefix.pdb.gz>Download receptor structure</a>]
    [<a href=?pdb=$pdbid&chain=$asym_id target=_blank>View receptor structure</a>]
    </td></tr>
    <tr><td><span title="Only residues with experimentally determined coordinates are included. Residues unobserved in the structure are excluded.">$seq_txt</span></td></tr>
    </table>
</div>
</td></tr>
'''.replace("$pdbid",pdbid
    ).replace("$asym_id",asym_id
    ).replace("$prefix",prefix
    ).replace("$bs",bs
    ).replace("$L",str(len(sequence))
    ).replace("$sequence",sequence
    ).replace("$seq_txt",seq_txt
    ))

    
    cmd="zcat %s/data/lig_all.tsv.gz|grep -P '^%s\\t%s\\t%s\\t'"%(
        rootdir,pdbid,asym_id,bs)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    line     =stdout.decode().splitlines()[0]
    items    =line.split('\t')
    lig3     =items[3]
    ligCha   =items[4]
    ligIdx   =items[5]
    resOrig  =items[6]
    resRenu  =items[7]
    manual   =items[8]
    moad     =items[9]
    pdbbind  =items[10]
    bindingdb=items[11]

    baff_list=[]
    if manual:
        baff_list.append("Manual survey: "+manual)
    if moad:
        baff_list.append("<a href=http://bindingmoad.org/pdbrecords/index/%s target=_blank>MOAD</a>: %s"%(pdbid,moad))
    if pdbbind:
        baff_list.append("<a href=http://pdbbind.org.cn/quickpdb.php?quickpdb=%s target=_blank>PDBbind-CN</a>: %s"%(pdbid,pdbbind))
    if bindingdb:
        baff_list.append("BindingDB: "+bindingdb)

    lig_prefix="%s_%s_%s_%s"%(pdbid,lig3,ligCha,ligIdx)
    filename="%s/output/%s.pdb.gz"%(rootdir,lig_prefix)
    x_list=[]
    y_list=[]
    z_list=[]
    if not os.path.isfile(filename):
        divided=pdbid[-3:-1]
        tar = tarfile.open("%s/weekly/ligand_%s.tar.bz2"%(rootdir,divided))
        fin=tar.extractfile("ligand/%s.pdb"%lig_prefix)
        lines=fin.read().decode().splitlines()
        txt=''
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" L"+line[22:]+'\n'
                x_list.append(float(line[30:38]))
                y_list.append(float(line[38:46]))
                z_list.append(float(line[46:54]))
            else:
                txt+=line+'\n'
        fout=gzip.open(filename,'wt')
        fout.write(txt)
        fout.close()
        fin.close()
        tar.close()
    else:
        fin=gzip.open(filename,'rt')
        for line in fin.read().splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x_list.append(float(line[30:38]))
                y_list.append(float(line[38:46]))
                z_list.append(float(line[46:54]))
        fin.close()
    xcen=sum(x_list)/len(x_list)
    ycen=sum(y_list)/len(y_list)
    zcen=sum(z_list)/len(z_list)
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Ligand information</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    ''')
    if lig3 in ["peptide","rna","dna"]:
        code=lig3.upper()
        if lig3=="peptide":
            code=lig3
        cmd="zcat %s/data/%s.fasta.gz|grep -PA1 '^>%s_%s_%s\\t'|tail -1"%(
            rootdir,lig3,pdbid,lig3,ligCha)
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        lig_sequence=stdout.decode().strip()
        lig_seq_txt='<br>'.join(textwrap.wrap(lig_sequence,50))
        print('''
<tr><td>&gt;$pdbid Chain $asym_id (length=$L)
[<a href=ssearch.cgi?seq_type=$lig3&sequence=$sequence>Search $code sequence</a>]
[<a href=output/$prefix.pdb.gz>Download ligand structure</a>]
[<a href=?pdb=$pdbid&chain=$asym_id&lig3=$lig3&ligIdx=0 target=_blank>View ligand structure</a>]
</td></tr>
<tr><td><span title="Only residues with experimentally determined coordinates are included. Residues unobserved in the structure are excluded.">$seq_txt</span></td></tr>
      '''.replace("$pdbid",pdbid
        ).replace("$asym_id",ligCha
        ).replace("$prefix",lig_prefix
        ).replace("$L",str(len(lig_sequence))
        ).replace("$lig3",lig3
        ).replace("$code",code
        ).replace("$sequence",lig_sequence
        ).replace("$seq_txt",lig_seq_txt
        ))
    else:
        formula=''
        InChI=''
        InChIKey=''
        SMILES=''
        name=''
        cmd="zcat %s/data/ligand.tsv.gz|grep -P '^%s\\t'"%(rootdir,lig3)
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        items=stdout.decode().split('\t')
        if len(items)>=6:
            formula,InChI,InChIKey,SMILES,name=items[1:6]
        svg="https://cdn.rcsb.org/images/ccd/labeled/%s/%s.svg"%(lig3[0],lig3)
        print('''
<tr>
    <td>
        <a href=$svg target=_blank><img src=$svg width=400></a>
    </td>
    <td>
    <table>
        <tr BGCOLOR="#DEDEDE"><td align=center><a href=https://wwpdb.org/data/ccd target=_blank>Ligand ID</a></td><td><span title="The ligand ID follows the
Chemical Component Dictionary (CCD)
used by the PDB database."><a href=https://rcsb.org/ligand/$lig3>$lig3</a></span></td></tr>
        <tr><td align=center><a href=https://inchi-trust.org target=_blank>InChI</a></td><td>$InChI</td></tr>
        <tr BGCOLOR="#DEDEDE"><td align=center><a href=https://inchi-trust.org target=_blank>InChIKey</a></td><td>$InChIKey</td></tr>
        <tr><td align=center>SMILES</td><td>$SMILES</td></tr>
        <tr BGCOLOR="#DEDEDE"><td align=center>Formula</td><td>$formula</td></tr>
        <tr><td align=center>Name</td><td>$name</td></tr>
        <tr BGCOLOR="#DEDEDE"><td align=center>Chain</td><td>$pdbid Chain $asym_id [<a href=output/$prefix.pdb.gz>Download ligand structure</a>] [<a href=?pdb=$pdbid&chain=$asym_id&lig3=$lig3&ligIdx=$ligIdx target=_blank>View ligand structure</a>]
        </td></tr>
    </table>
    </td>
</tr>
      '''.replace("$pdbid",pdbid
        ).replace("$asym_id",ligCha
        ).replace("$lig3",lig3
        ).replace("$ligIdx",ligIdx
        ).replace("$prefix",lig_prefix
        ).replace("$svg",svg
        ).replace("$InChIKey",InChIKey
        ).replace("$InChI",InChI
        ).replace("$SMILES",SMILES.replace(';',';<br>')
        ).replace("$formula",formula
        ).replace("$name",name.replace(';',';<br>')
        ))
    print('''
    </table>
</div>
</td></tr>''')
        
    reso   =''
    csaOrig=''
    csaRenu=''
    ec     =''
    go     =''
    uniprot=''
    pubmed =''
    cmd="zcat %s/data/pdb_all.tsv.gz |grep -P '^%s\\t%s\\t'|head -1"%(
        rootdir,pdbid,asym_id)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    items=stdout.decode().split('\t')
    if len(items)>=9:
        reso,csaOrig,csaRenu,ec,go,uniprot,pubmed=items[2:9]
        if reso=="-1.00":
            reso="N/A"
        else:
            reso=reso+" &#8491;"
    filename="%s/output/%s.pdb.gz"%(rootdir,prefix)
    if not os.path.isfile(filename):
        divided=pdbid[-3:-1]
        tar = tarfile.open("%s/weekly/receptor_%s.tar.bz2"%(rootdir,divided))
        fin=tar.extractfile("receptor/%s.pdb"%prefix)
        lines=fin.read().decode().splitlines()
        txt=''
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                txt+=line[:20]+" R"+line[22:]+'\n'
            else:
                txt+=line+'\n'
        fout=gzip.open(filename,'wt')
        fout.write(txt)
        fout.close()
        fin.close()
        tar.close()

    script=''
    if resOrig:
        resi_list=[r[1:] for r in resOrig.split()]
        script="select "+','.join(resi_list)+"; spacefill 25%; wireframe 50; color group;"
        for resi in resi_list:
            script+="select "+resi+" and *.ca; label %m%R; color label magenta;"
    script_ligand='spacefill 70% '
    if lig3=="peptide":
        script_ligand='cartoon; spacefill off;'
    if lig3 in ["rna","dna"]:
        script_ligand='cartoons; color grey; spacefill off; wireframe off'
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">Receptor-Ligand Complex Structure</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr BGCOLOR="#DEDEDE" align=center>
        <th>Global view</th><th>Local view</th><th>Structure summary</th>
    </tr>
    <tr>
        <td align=center>
<script type="text/javascript"> 
$(document).ready(function()
{
    Info = {
        width: 400,
        height: 400,
        j2sPath: "jsmol/j2s",
        script: "load output/$prefix.pdb.gz; color background black; cartoons; color group; spacefill off; wireframe off; $script; load append output/$lig_prefix.pdb.gz; select hetero; $script_ligand; frame all;"
    }
    $("#mydiv0").html(Jmol.getAppletHtml("jmolApplet0",Info))
});
</script>
<span id=mydiv0></span>
<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'spin on')">Spin on</a>]
[<a href="javascript:Jmol.script(jmolApplet0, 'spin off')">Spin off</a>]
[<a href="javascript:Jmol.script(jmolApplet0, 'Reset')">Reset</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay true')">High quality</a>]
[<a href="javascript:Jmol.script(jmolApplet0, 'set antialiasDisplay false')">Low quality</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet0, 'color background white')">White quality</a>]
[<a href="javascript:Jmol.script(jmolApplet0, 'color background black')">Black quality</a>]
        </td>
        <td align=center>
<script type="text/javascript"> 
$(document).ready(function()
{
    Info = {
        width: 400,
        height: 400,
        j2sPath: "jsmol/j2s",
        script: "load output/$prefix.pdb.gz; color background black; color group; spacefill off; wireframe off; $script; load append output/$lig_prefix.pdb.gz; select hetero; $script_ligand; zoomto 0 {$xcen $ycen $zcen}; frame all;"
    }
    $("#mydiv1").html(Jmol.getAppletHtml("jmolApplet1",Info))
});
</script>
<span id=mydiv1></span>
<br>
[<a href="javascript:Jmol.script(jmolApplet1, 'spin on')">Spin on</a>]
[<a href="javascript:Jmol.script(jmolApplet1, 'spin off')">Spin off</a>]
[<a href="javascript:Jmol.script(jmolApplet1, 'Reset')">Reset</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet1, 'set antialiasDisplay true')">High quality</a>]
[<a href="javascript:Jmol.script(jmolApplet1, 'set antialiasDisplay false')">Low quality</a>]<br>
[<a href="javascript:Jmol.script(jmolApplet1, 'color background white')">White quality</a>]
[<a href="javascript:Jmol.script(jmolApplet1, 'color background black')">Black quality</a>]
        </td>
        <td>
        <table>
            <tr><td align=center><strong>Resolution</strong><td>$reso</td></tr>
            <tr BGCOLOR="#DEDEDE"><td align=center><strong>Binding residue<br>(original residue number in PDB)</strong><td>$resOrig</td></tr>
            <tr><td align=center><strong>Binding residue<br>(residue number reindexed from 1)</strong><td>$resRenu</td></tr>
            <tr BGCOLOR="#DEDEDE"><td align=center><strong>Binding affinity</strong><td>$baff</td></tr>
        </table>
        </td>
    </tr>
    </table>
</div>
</td></tr>
    '''.replace("$pdbid",pdbid
    ).replace("$asym_id",asym_id
    ).replace("$reso",reso
    ).replace("$lig_prefix",lig_prefix
    ).replace("$prefix",prefix
    ).replace("$script_ligand",script_ligand
    ).replace("$script",script
    ).replace("$resOrig",resOrig
    ).replace("$resRenu",resRenu
    ).replace("$baff",'<br>'.join(baff_list)
    ).replace("$xcen","%.3f"%xcen,
    ).replace("$ycen","%.3f"%ycen,
    ).replace("$zcen","%.3f"%zcen,
    ))

    if ec:
        display_ec(ec,csaOrig,csaRenu)
    if go:
        display_go(go)
    return pubmed,uniprot

if __name__=="__main__":
    form   =cgi.FieldStorage()
    pdbid  =form.getfirst("pdb",'').lower()
    if not pdbid:
        pdbid=form.getfirst("pdbid",'').lower()
    asym_id=form.getfirst("chain",'')
    bs     =form.getfirst("bs",'')
    ligIdx =form.getfirst("idx",'')
    lig3   =form.getfirst("lig3",'')

    print("Content-type: text/html\n")
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/BioLiP1.png ></br>

<script type="text/javascript" src="jsmol/JSmol.min.js"></script>
<table style="table-layout:fixed;" width="100%" cellpadding="2" cellspacing="0">
<table>
''')
    pubmed=''
    uniprot=''
    if bs.startswith("BS"):
        pubmed,uniprot=display_interaction(pdbid,asym_id,bs)
    else:
        if lig3:
            if lig3 in ["rna","dna","peptide"]:
                pubmed=display_polymer_ligand(pdbid,asym_id,lig3)
            else:
                pubmed=display_regular_ligand(pdbid,asym_id,lig3,ligIdx)
        else:
            pubmed,uniprot=display_protein_receptor(pdbid,asym_id)
    
    uniprot_line=''
    if uniprot:
        uniprot_line='    <tr BGCOLOR="#DEDEDE"><td align=center><strong>UniProt</strong></td><td>'
        for u in uniprot.split(','):
            uniprot_line+="<a href=https://uniprot.org/uniprot/"+ \
                u+" target=_blank>"+u+"</a>, "
        uniprot_line=uniprot_line[:-2]+"</tr>\n"
        
    print('''
<tr><td>
<div id="headerDiv">
    <div id="titleText">External link</div>
</div>
<div style="clear:both;"></div>
<div id="contentDiv">
    <div id="RContent" style="display: block;">
    <table width=100% border="0" style="font-family:Monospace;font-size:14px;background:#F2F2F2;" >
    <tr><td align=center width=10%><strong>PDB</strong></td>
        <td width=90%><a href=https://rcsb.org/structure/$pdbid target=_blank>RCSB:$pdbid</a>,
            <a href=https://ebi.ac.uk/pdbe/entry/pdb/$pdbid target=_blank>PDBe:$pdbid</a>,
            <a href=https://pdbj.org/mine/summary/$pdbid target=_blank>PDBj:$pdbid</a></td> </tr>
    <tr BGCOLOR="#DEDEDE"> <td align=center><strong>PDBsum</strong></td><td><a href=http://ebi.ac.uk/pdbsum/$pdbid target=_blank>$pdbid</a></td></tr>
    <tr> <td align=center><strong>PubMed</strong></td><td><a href=https://pubmed.ncbi.nlm.nih.gov/$pubmed target=_blank>$pubmed</a></td></tr>
    $uniprot_line
    </tr>
'''.replace("$pdbid",pdbid
  ).replace("$pubmed",pubmed
  ).replace("$uniprot_line",uniprot_line
  ))
    
    print('''</table>
<p></p>
[<a href=.>Back to BioLiP</a>]
</body> </html>''')

#!/bin/bash

# This folder is the temporary folder for writing pdb files by cgi script.
# Permission and context of the folder and cgi scripts must be correctly set:

# First, set the permission:

    chmod a+x index.cgi
    chmod a+x ligand.cgi
    chmod a+x pdb.cgi
    chmod a+x qsearch.cgi
    chmod a+x ssearch.cgi
    chmod a+x esearch.cgi
    chmod a+x sym.cgi
    chmod a+x graphviz/bin/dot
    chmod 777 output/

# Second, set the context, which may not be necessary on some system:

    chcon -t httpd_sys_script_exec_t index.cgi
    chcon -t httpd_sys_script_exec_t ligand.cgi
    chcon -t httpd_sys_script_exec_t pdb.cgi
    chcon -t httpd_sys_script_exec_t qsearch.cgi
    chcon -t httpd_sys_script_exec_t ssearch.cgi
    chcon -t httpd_sys_script_exec_t esearch.cgi
    chcon -t httpd_sys_script_exec_t sym.cgi
    chcon -t httpd_sys_script_exec_t script/receptor1
    chcon -t httpd_sys_script_exec_t script/NWalign
    chcon -t httpd_sys_script_exec_t script/blastn
    chcon -t httpd_sys_script_exec_t script/blastp
    chcon -t httpd_sys_script_exec_t graphviz/bin/dot
    chcon -t httpd_sys_script_exec_t graphviz/lib/*so*
    chcon -t httpd_sys_script_exec_t graphviz/lib/graphviz/*so*
    chcon -t httpd_sys_rw_content_t  output/

# You can check the context by

    ls -laZ *cgi script/NWalign script/blastn script/blastp graphviz/bin/dot script/receptor1
    ls -laZ|grep output

#!/bin/bash
# skipde.sh LOG: per accepted cycle increment of the cumulative floor_audit [update] skip_de (last [update] line of each cycle = the accepted one; difference to the previous cycle's)
grep 'floor_audit cycle' "$1" | grep '\[update\]' | awk '{c=$4+0; for(i=1;i<=NF;i++) if($i=="skip_de"){v[c]=$(i+1)+0}; if(!(c in seen)){seen[c]=1; o[++n]=c}} END{p=""; for(k=1;k<=n;k++){c=o[k]; if(p!="") printf "%d %.6e\n", c, v[c]-v[p]; else printf "%d %.6e (first: cumulative value)\n", c, v[c]; p=c}}'

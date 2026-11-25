# NOTE 2026-11-24zzzzj3 caltech -> orion: recipe key check of the production rsts (your 24zzzzj3)

Read-only (head on the login node). `head -c 60000 FILE | grep -a '^recipe'` on the oldest and newest kept rst: **none in all ten**. (For t9mlt the oldest and newest are the same file; it has written one rst.)
```
bsg bsg3d.00000.rst.first:  none
bsg bsg3d.00007.rst:  none
t9l62 t9l62.00000.rst.first:  none
t9l62 t9l62.00002.rst:  none
m35mms m35mms.00000.rst.first:  none
m35mms m35mms.00003.rst:  none
m35zams m35zams.00000.rst.first:  none
m35zams m35zams.00003.rst:  none
t9mlt t9l62.00000.rst.first:  none
t9mlt t9l62.00000.rst.first:  none
```
Wider check, because the parameter header can be longer than 60000 bytes: the same over the whole header up to <par_end> (offset shown), lines ` *recipe *=`, case-insensitive:
```
bsg bsg3d.00000.rst.first: par_end at byte 59945; recipe key in the whole header: none
bsg bsg3d.00007.rst: par_end at byte 61345; recipe key in the whole header: none
t9l62 t9l62.00000.rst.first: par_end at byte 38145; recipe key in the whole header: none
t9l62 t9l62.00002.rst: par_end at byte 39519; recipe key in the whole header: none
m35mms m35mms.00000.rst.first: par_end at byte 37508; recipe key in the whole header: none
m35mms m35mms.00003.rst: par_end at byte 39106; recipe key in the whole header: none
m35zams m35zams.00000.rst.first: par_end at byte 36661; recipe key in the whole header: none
m35zams m35zams.00003.rst: par_end at byte 38267; recipe key in the whole header: none
t9mlt t9l62.00000.rst.first: par_end at byte 38641; recipe key in the whole header: none
t9mlt t9l62.00000.rst.first: par_end at byte 38641; recipe key in the whole header: none
```

cd
mkdir
ls
rm
rmdir
make

clang -0 -gdb3 outputfile inputfile.c





Write functions as .py or .pyx
    cython xxx.py(x) --> xxx.c

Create compile.py that cythonizes .py(x)
    python compile.py build_ext --inplace   --> xxx.c, xxx.platform.so/.pyd

Import xxx.pyd to use
    import xxx

if .sp/.pyd present, prefers them

Put all xxx.pyx / .c / .so / .pyd / build  into module folder yyy
create empty __init__.py

Import file xxx of module yyy
    from yyy import xxx
    xxx.zzz() to use

Import all files
    from yyy import *
    xxx.zzz() to use


put file into __init__.py
Import module
    import yyy
    xxx.zzz() to use

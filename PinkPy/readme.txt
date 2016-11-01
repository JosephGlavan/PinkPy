Copy contents to appropriate directory (where you keep the source files for your other Python packages):
    ... MANIFEST.in
    ... pinkGen.h
    ... pinkGen.cpp
    ... pinkpymodule.cpp
    ... setup.py
    
Make sure the path in "MANIFEST.in" points to the include directory for NumPy. The compiler needs to be able to find the NumPy API for C++.

Run "sudo python setup.py install"
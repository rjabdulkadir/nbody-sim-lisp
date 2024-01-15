This is an n-body simulation in Lisp derived from Python code by Philip Mocz. It uses nodgui (Lisp interface to Tk) and numcl (a Lisp implementation like numpy). The plan is to try more performant GUI and matrix operation libraries.

To run:
1. Install Tk
2. Install SBCL.
3. Install Quicklisp
4. On the command line:
5.   sbcl --dynamic-space-size 3076 --load nbody-0.1.lisp

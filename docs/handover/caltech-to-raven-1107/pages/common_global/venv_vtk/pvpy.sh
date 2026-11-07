#!/bin/bash
# venv vtk-osmesa python with a clean PYTHONPATH (module py-matplotlib numpy leaks otherwise)
unset PYTHONPATH PYTHONHOME
exec /resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python "$@"

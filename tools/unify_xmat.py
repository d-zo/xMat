#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
unify_xmat.py   v0.1
2023 Dominik Zobel
"""


# -------------------------------------------------------------------------------------------------
def Read_And_Custom_Parse(basepath, filename):
    import os

    content = ''
    with open(basepath + os.sep + filename, 'r') as infile:
        for singleline in infile:
            if (singleline.startswith('#addfile')):
                # Extract name of next file
                parts = singleline.split()
                nextfilename = parts[2] + '.f'
                content += Read_And_Custom_Parse(basepath=basepath, filename=nextfilename)
            elif (singleline.startswith('#adddirectory')):
                # Extract type and folder name
                parts = singleline.split()
                nextfilename = '_' + parts[1] + '.f'
                content += Read_And_Custom_Parse(basepath=basepath + os.sep + parts[2],
                    filename=nextfilename)
            else:
                content += singleline

    return content


# -------------------------------------------------------------------------------------------------
def main():
    content = Read_And_Custom_Parse(basepath='src', filename='xmat_structure.f')
    with open('xmat.f', 'w') as outfile:
        outfile.write(content)


main()

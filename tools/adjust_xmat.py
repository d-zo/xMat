#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
adjust_xmat.py   v0.2
2022-2023 Dominik Zobel
"""


import sys


# -------------------------------------------------------------------------------------------------
def Read_First_Line(filename):
    with open(filename, 'r') as infile:
        for textline in infile:
            firstline = textline.strip()
            break

    return firstline


# -------------------------------------------------------------------------------------------------
def Replace_Entries_In_Text(entries, text):
    mod_text = text
    for substring, replacement in entries.items():
        mod_text = mod_text.replace(substring, replacement)

    return mod_text


# -------------------------------------------------------------------------------------------------
def Replace_Entries_In_File(entries, filename):
    content = ''
    with open(filename, 'r') as infile:
        for textline in infile:
            content += Replace_Entries_In_Text(entries=entries, text=textline)

    with open(filename, 'w') as outfile:
        outfile.write(content)


# -------------------------------------------------------------------------------------------------
def Provide_Entries(version, representation):
    import datetime

    now = datetime.datetime.now()
    date = now.strftime('%Y-%m-%d')

    if (representation.upper() == 'TENS3333'):
        rep_mat = '3, 3'
        rep_tens = '3, 3, 3, 3'
        nelmat = '9'
    elif (representation.upper() == 'MAT99'):
        rep_mat = '9'
        rep_tens = '9, 9'
        nelmat = '9'
    elif (representation.upper() == 'MAT66'):
        rep_mat = '6'
        rep_tens = '6, 6'
        nelmat = '6'
    else:
        print('# Error: representation >' + representation + '< not supported')
        return None

    return {
        '--version--': version,
        '--cDatum--': date,
        '__matrix__': rep_mat,
        '__tensor__': rep_tens,
        '__nelmat__': nelmat
    }


# -------------------------------------------------------------------------------------------------
def main(args):
    if (len(args) != 2):
        print('# Error: File name expected')
        return

    version = '0.858'
    repline = Read_First_Line(filename='src/Math_Operations/_module.f')
    representation = repline.split()[2]

    entries = Provide_Entries(version=version, representation=representation)
    if (entries is None):
        print('# Error: Providing entries failed')
        return

    Replace_Entries_In_File(entries=entries, filename=args[1])


main(sys.argv)

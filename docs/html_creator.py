#!/usr/bin/python3


base_path = 'html'


def Write_File(filename, title, basetext, options, filelinks=dict()):
    content = '''<!DOCTYPE html>

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
   <meta http-equiv="content-type" content="text/html; charset=utf-8" />
   <title>Interactive help on using the xMat</title>
   <meta name="description" content="Interactive help on using the xMat" />
   <meta name="author" content="Dominik Zobel" />
   <meta name="keywords" lang="de" content="xMat, interactive help, constitutive models" />
   <meta name="keywords" lang="en" content="xMat, interaktive Hilfe, Stoffmodelle" />
   <meta name="viewport" content="width=device-width; initial-scale=1.0" />
   <style>
html, body {
   padding: 0;
   margin: 0;
   width: 100%;
   height: 100%;
   font-family: Helvetica,sans-serif;
}
* {
  -webkit-box-sizing: border-box;
  -moz-box-sizing: border-box;
  box-sizing: border-box;
}
body {
   display: flex;
   flex-flow: column nowrap;
}
.title {
   width: 100%;
   padding: 1em;
   font-size: 2em;
   font-weight: bold;
   color: #404040;
   background-color: #dddddd;
}
.text {
   padding: 2em 2em;
   background-color: #ffffff;
   line-height: 1.2em;
}
div.selection {
   background-color: #eeeeee;
   flex: 1 1 auto;
}
div.selection ul {
   margin: 0;
   list-style: none;
   display: flex;
   flex-flow: row wrap;
   padding: 1em;
}
div.selection ul li {
   margin: 1em;
   padding: 2em 0;
}
.pane {
   padding: 2em;
   border: 1px solid black;
   color: black;
   background-color: #88cc88;
   text-decoration: none;
}
.pane:hover {
   color: #eeeeee;
   background-color: #66aa66;
}
p {
   margin: 0;
   padding: 0;
}
pre {
   padding: 0.5em;
   margin: 0.4em 4em;
   background-color: #f2f2f2;
}
code {
   white-space: normal;
   background-color: #f2f2f2;
}
pre > code {
   white-space: pre-wrap;
   word-wrap: break-word;
   word-break: break-all;
   position: relative;
}
/*pre.numberSource code {
   counter-reset: source-line 0;
}
pre.numberSource code > span {
   position: relative;
   counter-increment: source-line;
}
pre.numberSource code > span > a:first-child::before {
   content: counter(source-line);
   position: relative;
   left: -3.2em;
   text-align: right;
   vertical-align: baseline;
   border: none;
   display: inline-block;
   -webkit-touch-callout:
   none; -webkit-user-select: none;
   -khtml-user-select: none;
   -moz-user-select: none;
   -ms-user-select: none;
   user-select: none;
   padding: 0 0.2em;
   width: 2em;
   background-color: #ffffff;
   color: #bfbfbf;
}*/
.note {
   padding: 2em;
   border: 1px solid black;
   color: black;
   background-color: #cccccc;
   text-decoration: none;
}
.note:hover {
   color: #eeeeee;
   background-color: #404040;
}
.twocolumns {
   columns: 2;
   -webkit-columns: 2;
   -moz-columns: 2;
}
.header_files {
   font-weight: bold;
}
ul.filelist, ul.dotless {
   list-style-type: none;
}
.filelist li {
   padding: 0.1em 0.4em;
}
.filelist li a {
   font-weight: bold;
   padding: 0.1em 0.4em;
   color: #888888;
   text-decoration: none;
}
.filelist li a:hover {
   background-color: #dddddd;
   color: #323232;
}
   </style>
</head>
<body>
<div class="title">''' + title + '''</div>
<div class="text">''' + basetext
    if (filelinks != {}):
        filelink_keylist = sorted(list(filelinks.keys()))
        content += '''<br /><br /><span class="header_files">Linked files</span><br />
<div class="twocolumns"><ul class="filelist">'''
        for key in filelink_keylist:
            if (filelinks[key] == ''):
                continue

            content += '   <li><a href="files/' + filelinks[key] + '"><svg xmlns="http://www.w3.org/2000/svg" width="1em" height="0.76em" viewBox="0 0 100 76"><path d="M 6,25 V 70 H 94 V 25" style="fill:none;stroke:#bbbbbb;stroke-width:12;stroke-linecap:round;stroke-linejoin:round" /><path d="M 58,3.0000002 V 22 H 74 L 50,50 26,22 H 42 V 3 Z" style="fill:#666666;stroke:#666666;stroke-width:6;stroke-linecap:round;stroke-linejoin:round" /></svg> ' + key + '</a></li>\n'

        content += '</ul>'

    content += '''</div></div>
<div class="selection">
<ul class="selection">
'''
    for optiontext, linknotes in options:
        link = linknotes[0]
        if (len(linknotes) > 1):
            classes = ' '.join(linknotes[1:])
        else:
            classes = 'pane'

        content += '   <li><a class="' + classes + '" href="' + link + '.html">' + optiontext + '</a></li>\n'

    content += '''</ul>
</div>
</body>
</html>
'''
    with open(filename, 'w') as outfile:
        outfile.write(content)


def Load_JSON(filename):
    import json

    with open(filename, 'r') as infile:
        content = json.load(infile)

    return content


def Base_Path():
    import os

    return base_path + os.sep


def Adjust_Option_List(basekey, options, refkeylist):
    import os

    nonexisting = []
    tempkeylist = sorted(list(options.keys()))
    processed_options = []
    for tempkey in tempkeylist:
        templink = options[tempkey][0]
        processed_options += [[tempkey, options[tempkey]]]
        if templink not in refkeylist:
            nonexisting += [templink]

    if (basekey != 'index'):
        processed_options += [['Return to beginning', ['index', 'note']]]

    for templink in nonexisting:
        print('Currently ' + templink + ' does not exist')
        refkeylist += [templink]
        linkname = Base_Path() + templink + '.html'
        title = 'This page is not finished and some information is still missing'
        basetext = 'Developer note: ' + templink
        tempoptions = [['Return to beginning', ['index', 'note']]]
        Write_File(filename=linkname, title=title, basetext=basetext, options=tempoptions)

    return refkeylist, processed_options


def main():
    import copy

    structure = Load_JSON(filename='help_selector.json')
    basekeylist = list(structure.keys())
    refkeylist = copy.deepcopy(basekeylist)
    for basekey in basekeylist:
        linkname = Base_Path() + basekey + '.html'
        num_successes = 0

        try:
            title = structure[basekey]['header']
            num_successes += 1
        except:
            title = ''

        try:
            basetext = ''.join(structure[basekey]['text'])
            num_successes += 1
        except:
            basetext = ''

        if (num_successes == 0):
            print('Warning: No content for ' + basekey + ' - skipping')
            continue

        # The following are considered optional and don't count against successes
        try:
            temp_options = structure[basekey]['choices']
        except:
            temp_options = dict()

        refkleylist, options = Adjust_Option_List(basekey=basekey,
            options=temp_options, refkeylist=refkeylist)

        try:
            filelinks = structure[basekey]['files']
        except:
            filelinks = dict()

        Write_File(filename=linkname, title=title, basetext=basetext, options=options, filelinks=filelinks)


main()

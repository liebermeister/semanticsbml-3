#!/usr/bin/env python

import os,string,copy

def shorten_names_and_numbers(inputstring):
    splits = []
    for single in inputstring.replace(" ","_").split("_"):
        if len(single) == 0: continue # no empty spaces
        if single[0] in string.letters: # if it is a string, take the first letter
            splits.append( single[0] )
        else: # we have a number
            if "e" in single:
                # scientific notation
                epoint = single.index("e")
                ppoint = single.index(".")
                newsingle = single[:ppoint+3]+single[epoint:]
                splits.append( newsingle ) if len(newsingle) < len(single) else splits.append( single )
            elif "." in single:
                # comma notation check for significant numbers
                ppoint = single.index(".")
                spoint = 0
                for c in single:
                    if c in "-+.0": spoint = spoint + 1
                    else: break
                ppoint = max(spoint,ppoint)
                splits.append( single[:ppoint+3] )
            else:
                splits.append(single)
    return "".join(splits)

class html:
    "docstring for html"
    def __init__(self):
        self.prebody = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\">\n<html>\n \
        <head>\n<title>TIde output</title>\n<style type=\"text/css\">\n \
        <!--\n \
        th{ background-color:#FFFF00 }\n \
        a{ color:#000000; text-decoration: none }\n \
        -->\n \
        </style>\n</head>\n<body>\n"
        self.postbody = "\n</body>\n</html>"

    def __create_submatrix__(self,matrix):
        newmatrix = []
        for l1 in matrix:
            newl = []
            for l2 in l1:
                if type(l2) == type([]):
                    l2 = self.__create_submatrix__(l2)
                    ht = htmltable("listbelow",indent=2)
                    ht.set_matrix(l2)
                    newl.append( ht )
                else:
                    newl.append( l2 )
            newmatrix.append(newl)
        return newmatrix
        
    def set_matrix(self,matrix):
        """Matrix can be a 2d or 4d matrix"""
        self.matrix = self.__create_submatrix__(matrix)
        self.htmlmatrix = htmltable(caption="",indent=1)
        self.htmlmatrix.set_matrix(self.matrix)

    def __create_subdic__(self,matrix):
        newmatrix = []
        for l1 in matrix:
            newl = []
            for l2 in l1:
                if type(l2) == type({}):
                    title = list(l2.keys())[0]
                    m = l2[title]
                    l2 = self.__create_subdic__(m)
                    ht = htmltable(title,indent=2)
                    ht.set_matrix(l2)
                    newl.append( ht )
                else:
                    newl.append( l2 )
            newmatrix.append(newl)
        return newmatrix
        
    def set_dic(self,dic):
        self.title = list(dic.keys())[0]
        self.matrix = self.__create_subdic__( dic[self.title] )
        self.htmlmatrix = htmltable(caption=self.title,indent=1)
        self.htmlmatrix.set_matrix(self.matrix)

    def __str__(self):
        return self.prebody+str(self.htmlmatrix)+self.postbody
    def tofile(self,filename):
        f = open(filename,"w")
        f.write( str(self) )
        f.close()

class htmltable:
    def __init__(self,caption="",indent=0):
        self.indentchar = " "
        self.indent = indent
        self.prelines = ["<table border=\"1\">","<caption>",self.indentchar+caption,"</caption>"]
        self.postlines = ["</table>"]
    def is_table(self):
        return True
    def set_matrix(self,matrix):
        self.matrix = matrix
        try:
            self.do_colouring = not matrix[0][0].is_table()
        except AttributeError:
            self.do_colouring = True
            # find maximum and minimum
            submatrix = [x[1:] for x in self.matrix[:-1]]
            l=[0.0]
            for x in submatrix: l=l+x
            while None in l: l.remove(None)
            self.maxentry = max(l)
            self.minentry = min(l)
    def __get_colour__(self,value):
        if value:
            fvalue = float(value)
        else:
            fvalue = 0.0
        if fvalue >= 0:
            positive = True
            relative = fvalue / float(self.maxentry) if self.maxentry != 0 else fvalue
        else:
            positive = False
            relative = fvalue / float(self.minentry) if self.minentry != 0 else fvalue
        r2hex = hex(int(256-255*relative))[2:].upper()
        while len(r2hex) < 2:
            r2hex = "0"+r2hex
        if positive:
            hexstring = r2hex+"FF"+r2hex
        else:
            hexstring = "FF"+r2hex+r2hex
        return hexstring
        
    def __generate_internal_lines__(self):
        lines = []
        lenmatrix = len(self.matrix)
        for row,index in zip(self.matrix,list(range(lenmatrix))):
            lines.append(self.indentchar+"<tr>")
            firstinrow = True
            for elem,index2 in zip(row,list(range(len(row)))):
                attributes = {}
                elemstring = str(elem) if elem != None else ""
                if self.do_colouring and (firstinrow or index == lenmatrix-1):
                    thistype = "th"
                    prehref = ""
                    posthref = ""
                    if index == lenmatrix-1:
                        # last row
                        elemstring = "<br/>".join([x for x in elemstring])
                else:
                    thistype = "td"
                    if self.do_colouring:
                        attributes["bgcolor"] = "#"+self.__get_colour__(elem)
                    href="/"+self.matrix[-1][index2]+"#"+self.matrix[index][0]+"="+elemstring
                    prehref = "<a href=\""+href+"\">"
                    posthref = "</a>"
                    if elemstring == "0.0" or elemstring == "0":
                        elemstring = ""
                # after the reference has been build, shorten the strings for the table
                #elemstring = shorten_names_and_numbers( elemstring )
                attributesstring = " ".join([ x+"=\""+y+"\"" for x,y in attributes.items()])
                lines.append(self.indentchar*2+"<"+thistype+" "+attributesstring+">"+prehref+elemstring+posthref+"</"+thistype+">")
                firstinrow = False
            lines.append(self.indentchar+"</tr>")
        return lines
    def __str__(self):
        internallines = self.__generate_internal_lines__()
        return "\n".join( [self.indentchar*self.indent+line for line in self.prelines+internallines+self.postlines] )


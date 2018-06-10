###################### 
#!/usr/bin/python
import numpy as np
import math
import re

################
filename='temp.top'
alpha = 0.4
################

def list_section(filename):
    fo = open(filename,'r')
    ls_section = []
    for line in fo:
        if re.match(r'^[[]',line):
             ls_section.append(line)
    return ls_section

def print_head(inputfile):
    fo = open(inputfile,'r')
    section=""
    for line in fo:
        if line.startswith('[ atomtypes ] ; 2'):
            break
        else:
            section += line
    fo.close()
    return section

def print_end(inputfile,start,end):
    fo = open(inputfile,'r')
    section=""
    flag='off'
    for line in fo:
        if line.startswith(start):
                #section += line
                flag='on'
        if flag=='on':
                section += line
        #if line.startswith(end):
        #    break
    fo.close()
    return section

def print_org(inputfile, start, stop):
    fo = open(inputfile,'r')
    section=""
    flag='off'
    for line in fo:
        if line.startswith(start):
            flag='on'
#            section += line
        if line.startswith(stop):
                #section += line
            flag='off'
            break;
        if flag=='on':
            section += line


    fo.close()
    return section

def read_section(inputfile, title, stop, scaling_factor, scaling_column):
    fo = open(inputfile,'r')
    section=""
    flag='off'
    for line in fo:
        if line.startswith(title):
            flag='on'
            section += line
            continue
        if flag=='on':
            if line.startswith(';'):
                section += line
                continue
            elif len(line.strip().split())==0:
                section += line
                #flag='off'
                #continue
            elif line.startswith(stop):
                #section += line
                flag='off'
                break;
            else:
                items = line.split()
                #print items
                scale = str(float(items[scaling_column])*scaling_factor)
                newitems = items[0:scaling_column] + [scale] + items[scaling_column+1:]
                newline = '\t'.join(newitems)+'\n'
                section += newline
    fo.close()
    return section

#def read_section_ildn(inputfile,title,stop,scaling_factor,scaling_column1, scaling_column2,scaling_column3):
def read_section_ildn(inputfile,title,stop,scaling_factor,scaling_column1, scaling_column2):
    fo = open(inputfile,'r')
    section=""
    flag='off'
    for line in fo:
        if line.startswith(title):
            flag='on'
            section += line
            continue
        if flag=='on':
            if line.startswith(';'):
                section += line
                continue
            elif len(line.strip().split())==0:
                section += line
                #flag='off'
                #continue
            elif line.startswith(stop):
                #section += line
                flag='off'
                break;
            else:
                items = line.split()
                #print items
                try:
                    scale1 = str(float(items[scaling_column1])*scaling_factor)
                    scale2 = str(float(items[scaling_column2])*scaling_factor)
                    #scale3 = str(float(items[scaling_column3])*scaling_factor)
                    #newitems = items[0:scaling_column1] + [scale1] + [scale2] + [scale3] + items[scaling_column1+3:]
                    newitems = items[0:scaling_column1] + [scale1] + [scale2] + items[scaling_column1+2:]
                    newline = '\t'.join(newitems)+'\n'
                    section += newline
                except IndexError:
                    section += line
    fo.close()
    return section

def read_section_lig(inputfile,title,stop,scaling_factor,scaling_column1, scaling_column2, scaling_column3, scaling_column4, scaling_column5, scaling_column6):
    fo = open(inputfile,'r')
    section=""
    flag='off'
    for line in fo:
        if line.startswith(title):
            flag='on'
            section += line
            continue
        if flag=='on':
            if line.startswith(';'):
                section += line
                continue
            elif len(line.strip().split())==0:
                section += line
                #flag='off'
                #continue
            elif line.startswith(stop):
                #section += line
                flag='off'
                break;
            else:
                items = line.split()
                #print items
                try:
                    scale1 = str(float(items[scaling_column1])*scaling_factor)
                    scale2 = str(float(items[scaling_column2])*scaling_factor)
                    scale3 = str(float(items[scaling_column3])*scaling_factor)
                    scale4 = str(float(items[scaling_column4])*scaling_factor)
                    scale5 = str(float(items[scaling_column5])*scaling_factor)
                    scale6 = str(float(items[scaling_column6])*scaling_factor)
                    #newitems = items[0:scaling_column1] + [scale1] + [scale2] + [scale3] + items[scaling_column1+3:]
                    newitems = items[0:scaling_column1] + [scale1] + [scale2] + [scale3] + [scale4] + [scale5] + [scale6] + items[scaling_column1+6:]
                    newline = '\t'.join(newitems)+'\n'
                    section += newline
                except IndexError:
                    section += line
    fo.close()
    return section

def read_section_ions(inputfile,title,stop,scaling_factor,scaling_column):
    fo = open(inputfile,'r')
    section=""
    flag='off'
    for line in fo:
        if line.startswith(title):
            flag='on'
            section += line
            continue
        if flag=='on':
            if line.startswith(';'):
                section += line
                continue
            elif len(line.strip().split())==0:
                section += line
                #flag='off'
                #continue
            elif line.startswith(stop):
                #section += line
                flag='off'
                break;
            else:
                items = line.split()
                #print items
                try:
                    scale = str(float(items[scaling_column])*scaling_factor)
                    newitems = items[0:scaling_column] + [scale] + items[scaling_column+1:]
                    newline = '\t'.join(newitems)+'\n'
                    section += newline
                except IndexError:
                    section += line
    fo.close()
    return section

def endline(inputfile):
    with open(inputfile, 'r') as f:
        return f.readlines()[-1]


sections = list_section(filename)
el = endline(filename)

print print_head(filename)
print read_section(filename, sections[1], sections[2],alpha,6)
print read_section(filename,sections[2],sections[3],alpha,4)
print read_section(filename,sections[3], sections[4], alpha,3)
print read_section(filename,sections[4], sections[5],alpha,5)
print read_section(filename,sections[5],sections[6],alpha,6)
print read_section(filename,sections[6],sections[7],alpha,6)
print print_org(filename,sections[7],sections[8])
print read_section(filename,sections[8],sections[9],alpha,6)

### chain A###
print print_org(filename,sections[9],sections[10])
print read_section(filename,sections[10], sections[11], math.sqrt(alpha),6)
print print_org(filename,sections[11],sections[14])
print read_section_ildn(filename, sections[14], sections[15], alpha, 5, 6)
print print_org(filename, sections[15],sections[16])
print print_org(filename, sections[16],sections[17])
#
##### chain B###
print print_org(filename,sections[17],sections[18])
print read_section(filename,sections[18], sections[19], math.sqrt(alpha),6)
print print_org(filename,sections[19],sections[22])
print read_section_ildn(filename, sections[22], sections[23], alpha, 5, 6)
print print_org(filename, sections[23],sections[25])
#
#### chain C###
print print_org(filename,sections[25],sections[26])
print read_section(filename,sections[26], sections[27], math.sqrt(alpha),6)
print print_org(filename,sections[27],sections[30])
print read_section_ildn(filename, sections[30], sections[31], alpha, 5, 6)
print print_org(filename, sections[31],sections[33])
#
##### chain D###
print print_org(filename,sections[33],sections[34])
print read_section(filename,sections[34], sections[35], math.sqrt(alpha),6)
print print_org(filename,sections[35],sections[38])
print read_section_ildn(filename, sections[38], sections[39], alpha, 5, 6)
print print_org(filename, sections[39],sections[41])
#
##### chain E###
print print_org(filename,sections[41],sections[42])
print read_section(filename,sections[42], sections[43], math.sqrt(alpha),6)
print print_org(filename,sections[43],sections[46])
print read_section_ildn(filename, sections[46], sections[47], alpha, 5, 6)
print print_org(filename, sections[47],sections[49])
#
##LIG##
print print_org(filename, sections[49],sections[50])
print read_section(filename,sections[50], sections[51], math.sqrt(alpha),6)
print read_section(filename,sections[51], sections[52], alpha,4)
print print_org(filename,sections[52],sections[53])
print read_section(filename,sections[53],sections[54],alpha,5)
print read_section_lig(filename, sections[54], sections[55], alpha, 5, 6, 7, 8, 9, 10)
print read_section(filename,sections[55], sections[56], alpha,6)

##Water and ions##
print print_org(filename,sections[56],sections[57])
print read_section(filename,sections[57], sections[58], math.sqrt(alpha),6)
print print_org(filename,sections[58],sections[60])
print read_section_ions(filename,sections[60],sections[-2],math.sqrt(alpha),6)
print print_org(filename,sections[-2],sections[-1])
print print_end(filename,sections[-1],el)

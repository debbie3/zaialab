# predicted_glycopeptides()
# Using output from expasy digestion (tryptic_peptides), glycosylation site from NetNGlyc (nglyc_sites), and possible glycoforms and their masses
# (potential_glycoforms)
# to determine theoretical masses of peptides including glycopeptides

import os

os.chdir("C:\Users\deborah\Documents\_Work and School\BU\Zaia Lab Rotation October 2015\projects\predicted_glycopeptides_HA_20151006")

def predicted_glycopeptides():

    peptides = open("tryptic_peptides.csv", 'r').readlines()
    peplist = []
    for line in peptides:
        line = line.strip() #removes extra newline character
        row = line.split(',') #splits mass, mc, sequence into a list
        peplist.append(row) #adds pep_item into peplist--nested lists

    nglyc = open("NGlyc_sites.csv", 'r').readlines()
    nglyclist = []
    for line in nglyc:
        line = line.strip()
        row = line.split(',')
        nglyclist.append(row) #add nglyc_item into nglyclist--nested lists

    glycoforms = open("potential_glycoforms.csv", 'r').readlines()
    glycoformslist = []
    for line in glycoforms:
        line = line.strip()
        row = line.split(',')
        glycoformslist.append(row)

    glycopeptides = open("predicted_glycopeptides.csv", 'w')
        
    # for each line in expasy output, search to see if nglyc_site matches with any part of the tryptic_peptide sequence
    for pep in peplist:
        position = pep[1].split('-')
        for glycsite in nglyclist:
            if glycsite[1] > position[0] and glycsite[1] < position[1] and (glycsite[5] == '++' or glycsite[5] == '+++'):
                for glycoform in glycoformslist:
                    mass = float(pep[0]) + float(glycoform[1])
                    print pep[3] + ', ' + str(mass) + ', ' +  'glycopeptide ' + glycoform[0] + '\n'
            
        print pep[3] + ', ' + pep[0] + '\n'
            
    glycopeptides.close()

predicted_glycopeptides()

# predicted_glycopeptides()
# Using output from expasy digestion (tryptic_peptides), glycosylation site from NetNGlyc (nglyc_sites), and possible glycoforms and their masses
# (potential_glycoforms)
# to determine theoretical masses of peptides including glycopeptides



# need to add argv stuff to have user input csv files
# write to file

import sys, itertools, os


cwd_path = os.getcwd()

program_name = sys.argv[0]

def predicted_glycopeptides():

    nglyc = open(cwd_path + '\\' + sys.argv[2], 'r').readlines()
    nglyclist = []
    for line in nglyc:
        line = line.strip() #removes extra newline character
        row = line.split(',') # splits each line into a list
        nglyclist.append(row) #nglyclist columns: stuff, position, sequence, stuff, stuff, +++'s

    glycoforms = open(cwd_path + '\\' + sys.argv[3], 'r').readlines()
    glycoformslist = [] #a list of names and masses of glycoforms
    glycomass = [] #a list of just the masses of glycoforms
    for line in glycoforms:
        line = line.strip()
        row = line.split(',')
        glycoformslist.append(row) #glycoformslist columns: glycan name, mass
        glycomass.append(float(row[1])) #glycomass column: mass

    glycopeptides = open("predicted_glycopeptides.csv", 'w')
        
    peptides = open(cwd_path + '\\' + sys.argv[1], 'r').readlines()
    for line in peptides:
        line = line.strip() #removes extra newline character
        pep = line.split(',') #splits mass, mc, sequence into a list
        position = pep[1].split('-')

        i = 0 # number of glycsites counter
        for glycsite in nglyclist:
        
            if glycsite[1] > position[0] and glycsite[1] < position[1] and (glycsite[5] == '++' or glycsite[5] == '+++'):
                i = i + 1
        
        if i == 0:
            print pep[3] + ', ' + pep[0] + '\n' #no glycsites

        else:
            permutationlist = itertools.product(glycomass, repeat = i)


            for permutation in permutationlist:
                mass = float(pep[0]) # initial mass of peptide
                for glycoform_mass in permutation:
                                
                    mass = mass + glycoform_mass
                                
                print pep[3] + ', ' + str(mass) + ',' + 'glycopeptide ' + '\n'
            
    glycopeptides.close()

predicted_glycopeptides()

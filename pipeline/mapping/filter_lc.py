#!/usr/bin/env python2

# Importation des modules
import sys, os

# Lecture des fichiers
try:
	f = open(sys.argv[1], "r")
	data = f.readlines()
	f.close()
except:
	print("Erreur chargement "+data+"\n")
	sys.exit(0) # Arret execution programme

try:
	f = open(sys.argv[2], "r")
	bad = f.readlines()
	f.close()
except:
	print("Erreur chargement"+bad+"\n")
	sys.exit(0) # Arret execution programme

# Creation dictionnaire
#bad_d = dict()
#for bad_read in bad :
#	bad_d[bad_read]=1

# Creation sous fichier
newfile=""
bad = map(lambda s: s.strip(), bad)
for line in data :
	id=line.split()[0]
	if (id not in bad) :
		newfile+=line
#	if (id not in bad_d) :
#		newfile+=line
		
f = open(sys.argv[3], "w")
f.write(newfile)
f.close()


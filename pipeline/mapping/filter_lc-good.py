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
	sys.exit(0)

try:
	f = open(sys.argv[2], "r")
	good = f.readlines()
	f.close()
except:
	print("Erreur chargement"+good+"\n")
	sys.exit(0)

# Creation sous fichier
newfile=""
good = map(lambda s: s.strip(), good)
for line in data :
	id=line.split()[0]
	if (id in good) :
		newfile+=line
		
f = open(sys.argv[3], "w")
f.write(newfile)
f.close()


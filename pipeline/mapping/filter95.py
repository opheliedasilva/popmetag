#!/usr/bin/env python2

#Importation des modules
import sys, os

#
encode=sys.argv[1]
complete=sys.argv[2]

try:
	f = open(encode, "r")
	fcode = f.readlines()
	f.close()
except:
	print("Erreur chargement fichier : verifiez existence fichier et relancer\n")
	sys.exit(0) # Arret execution programme

try:
	f = open(complete, "r")
	fall = f.readlines()
	f.close()
except:
	print("Erreur chargement fichier : verifiez existence fichier et relancer\n")
	sys.exit(0) # Arret execution programme

newfile=""
for i in range(len(fcode)) :
	code=fcode[i].split()[9]
	n=len(code)
	m=code.count('=')
	
	if (m/float(n) >= 0.95)	:
		newfile+=fall[i]

f = open(sys.argv[3], "w")
f.write(newfile)
f.close()


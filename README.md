# Calcul-Energie-Potentiel-de-l-eau
Programme Python qui calcule l’énergie potentielle totale de systèmes de molécules d’eau.
1- DESCRIPTION DU PROGRAMME:
	
	Programme Python qui calcule l’énergie potentielle totale de systèmes de molécules d’eau, dans le vide ou dans une boite de simulation. Pour cela, notre programme calcul les énergies des interactions liées (2 Bond et 1 Angle pour chaque molécules d’eau) ainsi que les énergies des interactions non liées lorsque le système contient au moins 2 molécules d’eau (3^N interactions de Van der Waals et 3^N interactions électrostatiques, où N représente le nombre de molécules d’eau)

	- Aprés avoir lancé le programme, une fenetre graphique s'affiche. Dans le menu déroulant il faut choisir le fichier pdb du système d'eau que vous souhaitez analyser puis appuyer sur "AFFICHER RÉSULTATS".  
	- Le calcul sera instantané pour de petits systèmes comme water.pdb ou water2.pdb mais peut prendre quelques secondes (10 à 20) pour un système plus grand comme water_box.pdb et affichera à la fin les résultats dans les champs appropriés à chaque Energie en kJ/mol et Kcal/mol.
	- Pour pouvoir visualiser l'emplacement des molécules d'eau (oxygène), l'utilisateur peut appuyer sur le bouton "Afficher Plot 3D" 


2- LES MODULES REQUIS POUR FAIRE TOURNER LE PROGRAMME:
 	Il faut avoir installé les modules suivants:
		 1- tkinter.
		 2- matplotlib.pyplot
		 3- mpl_toolkits.mplot3d

3- LIGNE DE COMMANDE POUR LANCER LE PROGRAMME: 
	 - En ayant activé miniconda : python script_youcef_thomas.py 
	 - Sans avoir activé miniconda : python3 script_youcef_thomas.py 


4- LES FICHIERS INPUT:
	Les fichiers d'entrées sont des fichiers .pdb contenant les coordonnées x,y,z pour les 3 atomes (O, H1 et H2) des molécules d'eau. 
	Si les molécules d'eau sont contenues dans une boite de simulation (ex: water_box.pdb), on retrouvera en plus une ligne débutant par "CRYST1" contenant les dimensions x,y,z de la boîte.

5- FICHIER OUTPUT:
 	 Les résultats seront aussi stockés dans un fichier de sortie, qui est par défaut nommé « sortie.txt », vous pouvez modifier et choisir le nom de ce fichier en complétant la case prévue à cet effet.

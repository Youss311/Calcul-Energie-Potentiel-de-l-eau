import math
from tkinter import *
from tkinter import ttk
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


#-----------Paramètres du champs de force ------------#

# Bond
l_eq=0.9572 # Dist à l'équilibre OH-H en Angstrom
K_bond=450.000 # Cste de raideur en kcal.mol-1.A-2

#Angle
teta_eq= 104.52*(math.pi/180) # angle converti en radian
K_angle= 55.000 #Cste de force en kcal.mol-1.rad-2

#LJ
#Profondeur du puit pour atomes ij (kcal.mol-1)
EpsHH=0.0460
EpsOH=0.0836
EpsOO=0.1521
#Distance à laquelle l'energie est minimale (A)
RminHH=0.4490
RminOH=1.9927
RminOO=3.5364

#Coulomb
EPS_r= 1 #Cste dielectrique pour le vide
#Charges --> unité e (charge électronique)
Q_O=-0.834 
Q_H=0.417
#Facteur de conversion e (kcal.A.mol-1.e-2)
f=332.0716

#Cutoff VdW et électrostatique
CUTOFF=14
#--------------Fonctions------------------#

#Lecture du fichier pdb
# molécules : liste de dictionnaire de dictionnaire
#Ex : [ {'OH2': {'x':..,'y':..,'z':..}, 'H1': {'x':..,'y':..,'z':..},'H2':{'x':.,'y':..,'z':..} } 
#       {'OH2': {'x':..,'y':..,'z':..}, 'H1': {'x':..,'y':..,'z':..},'H2':{'x':.,'y':..,'z':..} } ]

def isbox(fichier):
    box=False 
    Length_x=Length_y=Length_z=0
    with open (fichier, "r") as f_in:
        for ligne in f_in:
            if ligne.startswith("CRYST1"):
                box=True
                ligne=ligne.split()
                Length_x = int(float(ligne[1]))
                Length_y= int(float(ligne[2]))
                Length_z = int(float(ligne[3]))
    return box, Length_x, Length_y, Length_z


def read_pdb(fichier):
    liste_water=[]
    dico_atoms={}

    with open (fichier,"r") as f_in:
        for ligne in f_in:
            if ligne.startswith("ATOM") or ligne.startswith("HETATM"):
                dico_coord={  #Dictionnaire de coordonnees
                "x": float(ligne[30:38].strip()), 
                "y": float(ligne[38:46].strip()), 
                "z": float(ligne[46:54].strip()) }

                dico_atoms[ligne[12:16].strip()] = dico_coord #Dictionnaire de dictionnaire 
                #Dictionnaire avec les atomes et pour chaque atome : dictionnaire de coordonnees
                dico_coord={}

                if len(dico_atoms)==3: # Une fois que le dictionnaire contient les 3 atomes de la molécule d'eau 
                    # On ajoute dans une liste et on recommence la meme chose autant de fois qu'il y a de molécules
                    liste_water.append(dico_atoms)
                    dico_atoms={}

    return liste_water


# --- Champs de force --- #

# Fonction de liaison
def E_bond(l, l_eq, K_bond):
    return K_bond*((l-l_eq)**2)

#Fonction angle
def E_angle(teta, teta_eq, K_angle):
    return K_angle*((teta-teta_eq)**2)

# Fonction Lennard Jones (VdW)
def Lennard_Jones(Eps, Rmin, rij):
    return Eps*( (Rmin/rij)**12 -2*(Rmin/rij)**6 )

#Fonction Coulomb (Electrostatique)
def Coulomb(EPS_r, Qi, Qj, rij, f):
    return ((Qi*Qj)/(EPS_r*rij))*f


# --- Calcul termes lies --- #

#Calcul des vecteurs pour les angles
def calc_vec(OH, H):
    vec={
    "x": OH["x"]-H["x"],
    "y": OH["y"] - H["y"],
    "z": OH["z"] - H["z"]
    }
    return vec

#Calcul le produit scalaire de 2 vecteurs
def scalar(vecA,vecB):
    return (vecA['x']*vecB['x']) + (vecA['y']*vecB['y']) + (vecA['z']*vecB['z'])


#Calcul la norme de 2 vecteurs
def norme(vec):
    return math.sqrt(vec['x']**2 + vec['y']**2 +vec['z']**2 )


#Calcul des angles de valences
def angle(OH2,H1,H2):
    vecA=calc_vec(OH2,H1)
    vecB=calc_vec(OH2,H2)
    cos_teta=( scalar(vecA,vecB) / (norme(vecA) * norme(vecB)) ) 
    teta=math.acos(cos_teta)
    return teta


#Calcul de la distance entre 2 atomes (utile pour termes lies + non lies)
def dist(atmA,atmB):
    return math.sqrt( (atmA['x']-atmB['x'])**2 + (atmA['y']-atmB['y'])**2 + (atmA['z']-atmB['z'])**2 )


# Somme de tous les termes lies (angle + liaison)
def E_bonded_tot(syst_water):
    total_bonded = total_bond = total_angle = 0
    for water in syst_water: # On recupere les coordonnees (dictionnaire) de chaque atome
        for atoms in water:
            if atoms[0]=="O": #Si l'atome commence par O --> OH2 : dictionnaire contentant les coord de l'oxygene
                OH2=water[atoms]
            if atoms[0]=="H" and '1' in atoms: # Si l'atome commence par H et contient un 1 (Hydrogene1) --> H1 : dictionnaire contenant les coord de l'H1
                H1=water[atoms]
            else :  # Sinon (si pas O ou H1) --> c'est le H2 : dictionnaire contenant les coord de l'H2
                H2=water[atoms] 

        #Calcul des angles et de l'énergie
        teta=angle(OH2,H1,H2)
        Energie_angle=E_angle(teta,teta_eq,K_angle)
        total_angle += Energie_angle

        #Calcul des 2 termes de liaison et énergie
        l1=dist(H1,OH2)
        l2=dist(OH2,H2)
        Energie_bond= E_bond(l1, l_eq, K_bond) + E_bond (l2, l_eq, K_bond)

            #Somme des termes lies
        total_bond += Energie_bond

    return (total_angle, total_bond)


# --- Calcul termes non lies --- #


def PBC_dist (atmA, atmB, Length_x,Length_y,Length_z):
    half_x=Length_x/2
    half_y=Length_x/2
    half_z=Length_x/2

    xi,yi,zi = atmA.values()
    xj,yj,zj= atmB.values()

    dx = xi - xj#xj - xi

    if dx > half_x:
        dx -= Length_x
    if dx <= -half_x:
        dx += Length_x

    dy = yi - yj#yj - yi
    if dy > half_y:
        dy -= Length_y
    if dy <= -half_y:
        dy += Length_y
    dz = zi - zj#yj - yi
    if dz > half_z:
        dz -= Length_z
    if dz <= -half_z:
        dz += Length_z
    return math.sqrt( (dx)**2 + (dy)**2 + (dz)**2 )




def non_bonded_2waters(water1, water2, box, Length_x, Length_y, Length_z):
    LJ=Elec=r_00=0.0
    CUTOFF=14
    Calc_nobonded=False
    liste_r_HH=[] ; liste_r_OH=[]

    for atoms1 in water1: #On parcours toutes les interactions possibles entre atomes
        for atoms2 in water2:
            if box==True:  #Si présence d'une box : on prend en compte condition periodique aux limites pour calculer distance minimale
                rij=PBC_dist(water1[atoms1], water2[atoms2],Length_x,Length_y,Length_z)
            elif box == False: # Sinon on calcul juste la distance euclidienne normale
                rij=dist(water1[atoms1], water2[atoms2])

            if (atoms1[0]=="O" and atoms2[0]=="O"): # Si les 2 atomes sont des O --> on stocke la distance dans 1 variable (car 1 seule interaction O-O possible)
                r_OO=rij
            elif (atoms1[0]=="H" and atoms2[0]=="H"): # Si les 2 sont des H --> on stocke dans une liste (plusieurs interactions H-H possible (4))
                liste_r_HH.append(rij)
            else :
                liste_r_OH.append(rij) # Sinon : interaction O-H ou H-O --> on stocke dans une liste (plusieurs possibles (4))

            if box==True and rij<CUTOFF : #Pour la boite : Si une des interactions (sur les 9 possibles) est inferieur au cutoff
                Calc_nobonded=True # On met ce flag = True , et toutes les interactions non liés seront calculées
                #Si aucune des 9 distances n'est < 14 --> reste false et on calcul pas interaction non lié 

    if ((box==True and Calc_nobonded==True) or (box==False)) : #Cas ou on calcul les interaction 
        LJ+= Lennard_Jones(EpsOO, RminOO, r_OO ) #On calcul LJ et Elec pour l'interaction OO
        Elec+=Coulomb(EPS_r, Q_O, Q_O, r_OO, f)
        for r_HH in liste_r_HH : # On calcul LJ et Elec pour les interactions HH
            LJ+=Lennard_Jones(EpsHH, RminHH, r_HH)
            Elec+=Coulomb(EPS_r, Q_H, Q_H, r_HH, f)
        for r_OH in liste_r_OH: # On calcul LJ et Elec pour les interactions OH
            LJ+=Lennard_Jones(EpsOH, RminOH, r_OH)
            Elec+=Coulomb(EPS_r, Q_O, Q_H, r_OH, f)
        return LJ,Elec

    elif (box ==True and Calc_nobonded==False): # Si on a une boite est aucunes distance < 14 --> on return LJ=Elec=0
        return LJ, Elec # Dans ce cas la LJ et Elec valent tous les deux 0

# Somme de tous les termes non lies (VdW + Elec)
def sum_non_bonded(syst_water, box, Length_x,Length_y,Length_z):
    LJ_sum=Elec_sum=0
    for water1 in range(len(syst_water)-1): # Faire toutes les paires de molécules d'eau
        for water2 in range(water1+1, len(syst_water)):

            LJ,Elec=non_bonded_2waters(syst_water[water1],syst_water[water2],box,Length_x,Length_y,Length_z)
            LJ_sum+= LJ
            Elec_sum+=Elec

    return (LJ_sum, Elec_sum)

#============================================ Partie affichage ======================================================#
def PLOT_3D():
    fichier= entr_fich.get()
    fig=plt.figure()                       
    ax=fig.add_subplot(111,projection='3d')
    syst_water=read_pdb(fichier)           
    X=[] ; Y=[] ; Z=[]                     
    for water in syst_water: # On recupere les coordonnees (dictionnaire) de chaque atome
        for atoms in water:                
            if atoms[0]=="O": #Si l'atome commence par O --> OH2 : dictionnaire contentant les coord de l'oxygene
                X.append(water[atoms]["x"])
                Y.append(water[atoms]["y"])
                Z.append(water[atoms]["z"])
    ax.scatter(X,Y,Z, c='r',marker='o')
    ax.set_xlabel("X axis (Å)")
    ax.set_ylabel("Y axis (Å)")
    ax.set_zlabel("Z axis (Å)")
    ax.set_title("Plot 3D des positions des molécules d'eau (seulement les Oxygènes) ")
    plt.show()

# Fichier output
def fichier_output(nom_in,nom_fich_out,angl,bond,LJ,Elec,syst_water,Box,Length_x,Length_y,Length_z):
    elies=angl+bond
    enon_lies=LJ+Elec
    etot=angl+bond+Elec+LJ
    with open(nom_fich_out,"w") as f_out:
        f_out.write("Fichier analysé = {}\n\n".format(nom_in))  
        f_out.write("Nombre de molécules d'eau = {}\n".format(len(syst_water))) 
        f_out.write ("Résultats en Kcal/mol : \n" )
        f_out.write("E angle = {} Kcal/mol \n".format( angl))
        f_out.write("E bond = {} Kcal/mol \n".format (bond))
        f_out.write("E lies = {} Kcal/mol \n".format(elies))
        f_out.write("E VdW = {} Kcal/mol \n".format (LJ))
        f_out.write("E Elec = {} Kcal/mol \n".format (Elec))
        f_out.write("E non lies = {} Kj/mol \n".format (enon_lies))
        f_out.write("E total = {} Kcal/mol \n".format(etot))
        f_out.write("\n\n")
        f_out.write ("Résultats en Kj/mol : \n")
        f_out.write("E angle = {} Kj/mol \n".format( angl*4.184))
        f_out.write("E bond = {} Kj/mol \n".format (bond*4.184))
        f_out.write("E lies = {} Kj/mol \n".format((elies)*4.184))
        f_out.write("E VdW = {} Kj/mol \n".format (LJ*4.184))
        f_out.write("E Elec = {} Kj/mol \n".format (Elec*4.184))
        f_out.write("E non lies = {} Kj/mol \n".format ((enon_lies)*4.184))
        f_out.write("E total = {} Kj/mol \n".format((etot)*4.184))
        f_out.write("\n\n")
        if (Box==True):
            f_out.write("Une boite de simulation de {} Å x {} Å x {} Å est utilisée".format(Length_x,Length_y,Length_z) )
            

# Insertion des résultats dans les cases 
def Res():

    print("===================Calcul en cours=====================")

    nom_fich_out= "sortie.txt" # si l'utilisateur ne met pas de nom de fichier il aura sortie.txt comme nom standard 
    
    nom = entr_fich.get()
    syst_water=read_pdb(nom)
    Box, Length_x, Length_y, Length_z=isbox(nom)
    angl, bond= E_bonded_tot(syst_water)
    LJ, Elec = sum_non_bonded(syst_water, Box, Length_x,Length_y,Length_z)

    nbr_mol.delete(0,END) #supprime tous ce qui est dans une case avant un nouveau lancement 
    nbr_mol.insert(0,len(syst_water)) #insert le résultat dans une case 
    
    Angle.delete(0,END)
    Angle.insert(0,angl)
    Anglej.delete(0,END)
    Anglej.insert(0,angl*4.184)
    
    Bond.delete(0,END)
    Bond.insert(0,bond)
    Bondj.delete(0,END)
    Bondj.insert(0,bond*4.184)
    
    Lies.delete(0,END)
    Lies.insert(0,angl+bond)
    Liesj.delete(0,END)
    Liesj.insert(0,(angl+bond)*4.184)
    
    vdw.delete(0,END)
    vdw.insert(0,LJ)
    vdwj.delete(0,END)
    vdwj.insert(0,LJ*4.184)
    
    elec.delete(0,END)
    elec.insert(0,Elec)
    elecj.delete(0,END)
    elecj.insert(0,Elec*4.184)
    
    non_lies.delete(0,END)
    non_lies.insert(0,LJ+Elec)
    non_liesj.delete(0,END)
    non_liesj.insert(0,(LJ+Elec)*4.184)
    
    tot.delete(0,END)
    tot.insert(0,angl+bond+Elec+LJ)
    totj.delete(0,END)
    totj.insert(0,(angl+bond+Elec+LJ)*4.184)
    if (Box==True):
        box_in.delete(0,END)
        box_in.insert(0,"OUI")
        box_dim.delete(0,END)                                     
        box_dim.insert(0,"{} Å x {} Å x {} Å".format(Length_x , Length_y , Length_z))
    else:
        box_in.delete(0,END)
        box_in.insert(0,"NON")
        box_dim.delete(0,END)
        box_dim.insert(0,"XXXXXXXXXXXXXXXXXXX")


    print("===================Fin De Calcul=====================")

    if(fich_out.get()):
        nom_fich_out=fich_out.get()
    else:
        print("Vous n'avez pas mis de nom pour le fichier de sortie, ce dernier sera tout de même dans votre répertoire sous le nom : sortie.txt")
    fichier_output(nom,nom_fich_out,angl,bond,LJ,Elec,syst_water,Box,Length_x,Length_y,Length_z)

# Affichage de la fenetre 

def affichage(): 
    global entr_fich,Angle,Bond,Lies,vdw,elec,non_lies,tot,nbr_mol,Anglej,fich_out #ces variables (sont les cases 'entry') elle sont en global pour 
    global Bondj,Liesj,vdwj,elecj,non_liesj,totj,nbr_molj,box_in,box_dim  # que la fonction Res puisse les utiliser car on ne peux pas faire un return 
    #fenetre
    window=Tk()
    window.title("Calcul energie de molécules d'eau")
    window.geometry("700x650")
    window.minsize(800,650)
    window.maxsize(800,650)
    window.config(background='#859dd8')
    label_title=Label(window,fg='white',text="Bienvenue dans\n AQUA ENERGY CALCULATOR",font=("Helvetica",25),bg='#859dd8')
    label_title.pack()
    
    lbl_fich_pdb=Label(window,text="Choisissez le fichier pdb à analyser : ", font=("Helvetica",15),bg='#859dd8')
    lbl_fich_pdb.place(x=30,y=150)
    OptionList = ["water.pdb","water2.pdb","water3.pdb","water10.pdb","water_box.pdb","water3_bis.pdb","water_bis.pdb"]
    entr_fich = ttk.Combobox( window,values=OptionList)
    entr_fich.current(3) # pour afficher le premier element de la liste dans la combobox par défaut
    entr_fich.config(width=15, font=('Helvetica', 12))
    entr_fich.place(x=450,y=150)

    lbl_fich_out=Label(window,text="Entrez le nom souhaité pour le fichier de sortie : ", font=("Helvetica",15),bg='#859dd8')
    lbl_fich_out.place(x=30,y=185)
    fich_out=Entry(window,width=20)
    fich_out.place(x=460,y=190)
    
    lblbox=Label(window, text="Système avec BOX ? " ,bg='#859dd8',font=("Helvetica",15))
    lblbox.place(x=30,y=335)
    box_in=Entry(window,width=5)
    box_in.place(x=250,y=340)
    lbl_dim_box=Label(window, text="Dimensions de la BOX" ,bg='#859dd8',font=("Helvetica",15))
    lbl_dim_box.place(x=30,y=365)
    box_dim=Entry(window,width=20)
    box_dim.place(x=250,y=370)

    lbl_nbr_mol=Label(window, text="Le nombre de molécules d'eau = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_nbr_mol.place(x=30,y=395)
    nbr_mol=Entry(window)
    nbr_mol.place(x=330,y=400)
    lbl_mol=Label(window, text="Molécule(s)",bg='#859dd8',font=("Helvetica",15))
    lbl_mol.place(x=500,y=395)

    lbl_Angle=Label(window, text="E Angle = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_Angle.place(x=30,y=425)
    Angle=Entry(window,width=23)
    Angle.place(x=200,y=430)
    lbl1=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl1.place(x=400,y=425)
    Anglej=Entry(window,width=23)
    Anglej.place(x=500,y=430)
    lbl1j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl1j.place(x=700,y=425)

    lbl_Bond=Label(window, text="E Bond = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_Bond.place(x=30,y=455)
    Bond=Entry(window,width=23)
    Bond.place(x=200,y=460)
    lbl2=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl2.place(x=400,y=455)
    Bondj=Entry(window,width=23)
    Bondj.place(x=500,y=460)
    lbl2j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl2j.place(x=700,y=455)
    
    lbl_Lies=Label(window, text="E Liés = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_Lies.place(x=30,y=485)
    Lies=Entry(window,width=23)
    Lies.place(x=200,y=490)
    lbl3=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl3.place(x=400,y=485)
    Liesj=Entry(window,width=23)
    Liesj.place(x=500,y=490)
    lbl3j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl3j.place(x=700,y=485)
    
    lbl_vdw=Label(window, text="E Van der waals = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_vdw.place(x=30,y=515)
    vdw=Entry(window,width=23)
    vdw.place(x=200,y=520)
    lbl4=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl4.place(x=400,y=515)
    vdwj=Entry(window,width=23)
    vdwj.place(x=500,y=520)
    lbl4j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl4j.place(x=700,y=515)
    
    lbl_elec=Label(window, text="E Elec = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_elec.place(x=30,y=545)
    elec=Entry(window,width=23)
    elec.place(x=200,y=550)
    lbl5=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl5.place(x=400,y=545)
    elecj=Entry(window,width=23)
    elecj.place(x=500,y=550)
    lbl5j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl5j.place(x=700,y=545)
    
    lbl_non_lies=Label(window, text="E Non Liés = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_non_lies.place(x=30,y=575)
    non_lies=Entry(window,width=23)
    non_lies.place(x=200,y=580)
    lbl6=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl6.place(x=400,y=575)
    non_liesj=Entry(window,width=23)
    non_liesj.place(x=500,y=580)
    lbl6j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl6j.place(x=700,y=575)
    
    lbl_tot=Label(window, text="E Total = " ,bg='#859dd8',font=("Helvetica",15))
    lbl_tot.place(x=30,y=605)
    tot=Entry(window,width=23)
    tot.place(x=200,y=610)
    lbl7=Label(window, text="Kcal/mol =",bg='#859dd8',font=("Helvetica",15))
    lbl7.place(x=400,y=605)
    totj=Entry(window,width=23)
    totj.place(x=500,y=610)
    lbl7j=Label(window, text="Kj/mol",bg='#859dd8',font=("Helvetica",15))
    lbl7j.place(x=700,y=605)

    Result=Button(window, text="AFFICHER RESULTATS" ,font=('Helvetica', 25),fg='white' ,bg='#e08e06',command=Res)
    Result.place(x=100,y=250)
    plot= Button(window, text="Afficher Plot 3D" ,font=('Helvetica', 10),fg='white' ,bg='#e08e06',command=PLOT_3D)
    plot.place(x=550,y=265)
    
    window.mainloop()
    

# ===================================Main==========================================================================#

# -----------Main------------ #

if __name__ == "__main__":
    
    affichage()


    
#--------------------------------------------------------------




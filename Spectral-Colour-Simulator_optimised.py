import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize as sc
    
# Travail preliminaire
## Fichiers

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

file_path1="lampe.D65"
file_path2="D149_exp"
file_path3="oeil"
file_path4="lampe.A.txt"

print("------------------------------------------------------------------------ \n")
print("  SPECTRAL COLOUR SIMULATOR - 2025 \n")
print("------------------------------------------------------------------------ \n")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⡇⠀⠈⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣆⠀⠀⠀⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣷⣄⠀⠀⠀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣧⡀⠀⠀⠀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣦⡀⠀⠀⠀⠈⠙⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡄⠀⠀⠀⠀⠀⠈⣻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⡀⠀⢀⣴⠟⠉⠉⠻⣿⣿⣿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣾⠃⠀⠀⣠⠶⠛⠻⢿⣿⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⣀⣾⠁⠀⠀⠀⠀⠙⣿⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⠀⠀⠀⠀⠀⠀⢹⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣆⠀⠀⠀⠀⠀⢸⣿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⡀⠀⠈⠻⢿⣿⣿")
print("⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿ \n")
print("- Avant de commencer, assurez vous d'avoir mis tout les fichiers ('oeil', 'lampe.D65', 'lampe.A.txt')  dans le meme dossier que ce fichier .py")
print("- Si vous remplacez les donnees de ces fichiers (e.g. la source) le nombre de points experimentaux doit rester 400 \n")


def get_bool(prompt):
    while True:
        try:
            return {"oui": True, "non": False}[input(prompt).lower()]
        except KeyError:
            print ("Repondez oui ou non")
            
def get_int(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Tapez un nombre, la gamme recommandé se situe entre 0.1 et 2.5")

if not get_bool("Etes vous sur windows? \n"):
    file_path1="lampe.D65"
    file_path2="D149_exp"
    file_path3="oeil"
    file_path4="lampe.A.txt"
else:
    file_path1=os.path.join(__location__, "lampe.D65.txt")
    file_path2=os.path.join(__location__, "D149_exp.txt")
    file_path3=os.path.join(__location__, "oeil.txt")
    file_path4=os.path.join(__location__, "lampe.A.txt")

def main():
    while True:   
        if not get_bool("Utiliser les parametres de base (lumiere naturelle, sigma=0.5) ? - entrez oui ou non \n"):
            artificielle = get_bool("Vous voulez simuler une lumière artificielle - entrez oui ou non \n")
            sigma_global=get_int("Entrez la valeur de sigma souhaitée (valeur recommandé comprise entre 0.1 et 2.5) \n")
        else:
            artificielle=False
            sigma_global=0.5

        print("\n ------------------------------------------------------------------------ \n")

        ## Creation des listes de donnees
        Lambda_0=[]
        T_D149=[]
        S_D65=[]
        S_A=[]
        x_bar=[]
        y_bar=[]
        z_bar=[]


        with open(file_path1) as file:
            for _ in range(1):
                next(file)
            for line in file:
                a, b = line.rstrip().split("	", 1)
                Lambda_0.append(float(a))
                S_D65.append(float(b))

        with open(file_path2) as file:
            for _ in range(1):
                next(file)
            for line in file:
                a, b = line.rstrip().split(" ", 1)
                #Lambda2.append(float(a))
                T_D149.append(float(b)/100)

        with open(file_path3) as file:
            for _ in range(1):
                next(file)
            for line in file:
                a, b, c, d = line.rstrip().split("	")
                #Lambda.append(float(a))
                x_bar.append(float(b))
                y_bar.append(float(c))
                z_bar.append(float(d))

        with open(file_path4) as file:
            for _ in range(1):
                next(file)
            for line in file:
                a, b = line.rstrip().split("	", 1)
                S_A.append(float(b))


        T_D149.pop()
        x_bar.pop()
        y_bar.pop()
        z_bar.pop()


        ## Constantes
        lambda1=380 ##debut domaine spectral
        lambda2=779 ##fin du domaine

        C1=3.748406*10**(-16)
        C2=1.438769*10**(-2)
        C3=2885*10**(-6)
        h=6.62607015*10**(-34)
        c=299792458
        eV_J=1.6021774232052327*10**(-19)
        hc_eV = (h*c)/eV_J
        Lambda2=[i for i in range(lambda1, lambda2+1)]
        Lambda = [hc_eV/(x*1e-9) for x in Lambda2]
        ## Le visible est entre 1.55 et 3.1 eV
        print(Lambda)


        ## Creation des listes de donnees
        Excited = [1,2,3,4,5,6,7,8,9,10,11,12]
        Lambda_excite=[511.70, 473.63, 410.68, 399.84, 360.93, 357.87, 344.05, 331.97, 317.57, 313.58, 312.30, 301.93]
        f=[1.0296, 0.7305, 0.0291, 0.0022, 0.0834, 0.4263, 0.2625, 0.0200, 0.0751, 0.0013, 0.1172, 0.0155]
        E = [2.4230, 2.6177,3.0190,3.1008,3.4351,3.4645,3.6037,3.7348,3.9041,3.9538,3.9700,4.1064]

        ## Optimisation : Intervalle lambda initialision
        Excited2 = []
        Lambda_excite2 = []
        f2 = []
        E2 = []
        
        for i in range (len(Excited)):
            if lambda2>Lambda_excite[i]>lambda1-30 :
                print("ui")
                Excited2.append(Excited[i])
                Lambda_excite2.append(Lambda_excite[i])
                f2.append(f[i])
                E2.append(E[i])

        sigma = [sigma_global for i in range(len(E2))]

        ## Fonction Gaussiene

        def gaussienne(lambdax,lambdai,sigmai,f) :
            #return 2*np.sqrt(np.log(2)/np.pi)*(f/sigmai)*np.exp(-(2*(lambdax*lambdai-h*c)*np.sqrt(np.log(2))/(lambdai*sigmai))**2)
            return 2*np.sqrt(np.log(2)/np.pi)*(f/sigmai)*np.exp(-4*np.log(2)*((lambdax-lambdai)/sigmai)**2)

        def calcule_gi(lambdax,lambdai,sigmai,f,Scale) :
            Sum = 0
            if type(sigmai) is float:
                for i in range(len(lambdai)):
                    Sum+=gaussienne(lambdax,lambdai[i],sigmai,f[i])
            elif type(sigmai) is int:
                for i in range(len(lambdai)):
                    Sum+=gaussienne(lambdax,lambdai[i],sigmai,f[i])
            
            else :
                for i in range(len(lambdai)) :
                    Sum+=gaussienne(lambdax,lambdai[i],sigmai[i],f[i])
            Sum*=Scale
            return Sum

        ## On considere que f est une absorbance

        def transmittance(f):
            return 10**(-f)

        ## Residu

        def residuals_sigmas(params, Lambda, E2, f2, T_exp):
            # params = liste de sigma à optimiser (un par raie)
            G = [calcule_gi(x, E2, params, f2, 1) for x in Lambda]
            G_trans = np.array([transmittance(val) for val in G])
            return G_trans - (T_exp / max(T_exp))

        ## Gaussienne fitting
        tolerance = 0.01
        err = float('inf')
        max_iter = 15
        iteration = 0

        while err > tolerance and iteration < max_iter:
            result = sc.least_squares(
                residuals_sigmas,
                sigma,
                bounds=(0, np.inf),
                args=(Lambda, E2, f2, np.array(T_D149))
            )
            sigma = list(result.x)
            G = [calcule_gi(x, E2, sigma, f2, 1) for x in Lambda]
            G_trans = np.array([transmittance(val) for val in G])
            err = np.sqrt(np.mean((G_trans - np.array(T_D149)/max(T_D149))**2))
            iteration += 1
            print(f"Iteration {iteration}, erreur RMS = {err}")


        ## Gaussienne fitting
        G_non_opt=[]
        for x in Lambda:
            G_non_opt.append(calcule_gi(x,E,sigma_global,f,1))



        ### Graphes en eV et en lambda de l'absorbance
        #### Dirac
        #plt.scatter(Lambda_excite,f,s=3,c='blue')
        for i in range (len(E)): ## barres de type dirac
            plt.plot([E[i],E[i]],[f[i],0],color="blue")
        ### Gaus
        fig, ax = plt.subplots()
        plt.plot(Lambda, G, label="optimisé")
        plt.plot(Lambda, G_non_opt, label="non optimisé")
        plt.title("Spectre de Transmittance D149 theorique")
        plt.xlabel("Energie (eV)")
        plt.ylabel('Absorbance')
        ax.legend()
        plt.show()


        for i in range (len(E)): ## barres de type dirac
            plt.plot([Lambda_excite[i],Lambda_excite[i]],[f[i],0],color="blue")
            

        plt.plot(Lambda2, G)
        plt.title("Spectre Dirac et Gaus")
        plt.xlabel("Longeur d'onde (nm)")
        plt.ylabel('Absorbance')
        plt.show()

        # Graphe transmittance theorique et reel
        ## Theorie
        fig, ax = plt.subplots()

        a=max(T_D149)
        Trans=[transmittance(x)*a for x in G]
        plt.plot(Lambda2, Trans, label='T_theorique')
        plt.title("Spectres de transmission")
        plt.xlabel("Longeur d'onde (nm)")
        plt.ylabel('Transmittance')

        ## Reel
        plt.plot(Lambda_0, T_D149, label='T_reel')
        ax.legend()
        plt.show()



        # Couleurs reel et theorique
        ## Theorique
        ### Transmittance
        L_ech=[]
        for x in G:
            L_ech.append(transmittance(x))


        ### Fonction integration liste elements espaces de 1
        def integrale_rec_points(f): #methode rec de "fonction" (liste) sur son domaine
            aire=0
            for i in range (0,len(f)-1):
                aire+=(f[i]+f[i+1])/2
            return aire

        ### Creation listes integrandes theoriques
        integrande_x=[]
        integrande_y=[]
        integrande_z=[]
        integrande_K=[]

        def integrande(L):
            if artificielle==False:
                for i in range(lambda2-lambda1+1):
                   integrande_x.append(S_D65[i]*x_bar[i]*(L[i]))
                   integrande_y.append(S_D65[i]*y_bar[i]*(L[i]))
                   integrande_z.append(S_D65[i]*z_bar[i]*(L[i]))
                   integrande_K.append(S_D65[i]*y_bar[i])
            else:
                for i in range(lambda2-lambda1+1):
                   integrande_x.append(S_A[i]*x_bar[i]*(L[i]))
                   integrande_y.append(S_A[i]*y_bar[i]*(L[i]))
                   integrande_z.append(S_A[i]*z_bar[i]*(L[i]))
                   integrande_K.append(S_A[i]*y_bar[i])

        integrande(L_ech)


        ### Calcul de K, X, Y, Z par integration
        def integrale_para():
            K=100.0/integrale_rec_points(integrande_K)
            X=K*integrale_rec_points(integrande_x)
            Y=K*integrale_rec_points(integrande_y)
            Z=K*integrale_rec_points(integrande_z)
            return K, X, Y, Z

        K,X,Y,Z = integrale_para()
        print("X,Y,Z theoriques = ", [X, Y, Z])


        ### Changement base
        def changement_base(X, Y, Z):
            var_X = X / 100
            var_Y = Y / 100
            var_Z = Z / 100

            var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986
            var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415
            var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570

            if var_R > 0.0031308:
                var_R = 1.055 * ( var_R ** ( 1 / 2.4 ) ) - 0.055
            else:
                var_R = 12.92 * var_R
            if var_G > 0.0031308 :
                var_G = 1.055 * ( var_G ** ( 1 / 2.4 ) ) - 0.055
            else:
                var_G = 12.92 * var_G
            if var_B > 0.0031308:
                var_B = 1.055 * ( var_B ** ( 1 / 2.4 ) ) - 0.055
            else:
                var_B = 12.92 * var_B
            return var_R * 255, var_G * 255, var_B * 255

        sR, sG, sB=changement_base(X, Y, Z)

        ### Coupe les valeurs negatives

        def format_RBG(sR, sG, sB):
            sR=np.clip(sR, 0, 255)
            sG=np.clip(sG, 0, 255)
            sB=np.clip(sB, 0, 255)
            return sR, sG, sB
                   
        sR, sG, sB=format_RBG(sR, sG, sB)


        ### Affichage en cercle de couleur (theorique)
        fig = plt.figure()
        ax = fig.add_subplot()

        print("R, G, B theoriques =", sR, sG, sB)
        
        circ_T = plt.Circle((2,0), 4, color=(sR/255,sG/255,sB/255))
        ax.add_patch(circ_T)
        plt.text(1, 5, 'Theorique', fontsize=12, bbox=dict(facecolor='black', alpha=0))
        

        ### Calculs listes integrandes reels
        integrande_x=[]
        integrande_y=[]
        integrande_z=[]
        integrande_K=[]

        integrande(T_D149)

        K,X,Y,Z = integrale_para()
        print("X,Y,Z reels = ", [X, Y, Z])
        print("R, G, B reels =", sR, sG, sB)

        sR, sG, sB=changement_base(X, Y, Z)
        sR, sG, sB=format_RBG(sR, sG, sB)

        circ_R = plt.Circle((12,0), 4, color=(sR/255,sG/255,sB/255))
        ax.add_patch(circ_R)
        plt.text(11, 5, 'Réel', fontsize=12, bbox=dict(facecolor='black', alpha=0))


        ## ### Affichage en cercle de couleur (reel)

        plt.axis('equal')
        plt.axis('off')
        plt.show()

        if get_bool("Recommencer? \n"):
            pass
        else:
            print("------------------------------------------------------------------------ \n")
            print("Merci d'avoir utilisé SPECTRAL COLOUR SIMULATOR - 2025 \n")
            print("------------------------------------------------------------------------ \n")
            quit()

main()

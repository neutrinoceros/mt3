Résumé du travail accompli :
============================

#Analyse des données

- Données : MHB et REN dX,dY,erreurs
- interpolation (régularisation du signal)
- FFT : comparaison avec REN --> le pic à 450j est d'origine géophysique (evidence existence noyau) donc on doit l'éliminer
- obtention des amplitudes complexes et des phases des composantes (quelles composantes ??) du signal par moindre carrés
- rafinement : prise en compte d'un biais et d'une pente

#Traitement

- Elimination du terme géophysique une fois son amplitude connue puis reconstruction du signal nettoyé par FFT inverse
- rerun de l'ajustement par moindre carrés des amplitudes complexes
- ajustement de la fonction de transfert T(sig) = eta(sig)MHB / eta(sig)REN --> obtention de corrections aux paramètres géophysiques donnés a priori (Q : d'où sortaient les valeurs données ?)

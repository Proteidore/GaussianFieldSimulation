# GaussianFieldSimulation

English: This project is a set of benchmarks testing several methods of simulation of Gaussian fields.

Français: Ce projet est un ensemble de benchmarks testant plusieurs méthodes de simulation de champs gaussiens.

## Description de l'arborescence de fichiers

### Benchmarks

     ---> CholeskyMethod
     
          |__ cholesky.py   -> fait le benchmark de la méthode de Cholesky
          
          |__ testTool.py   ->  fonction pour le calcul de l'erreur L2
          
     ---> GalliGaoGibbsMethod
     
          |__ GGG.py   -> fait le benchmark de la méthode de Galli-Gao-Gibbs  
          
          |__ testTool.py   -> fonction pour le calcul de l'erreur L2
          
     ---> HmatrixMethod
     
          |__ Hmatrix.py   -> fait le benchmark de la méthode des H-matrices
          
          |__ testTool.py   -> fonction pour le calcul de l'erreur L2
          
     ---> P1InterpolationMethod
     
          |__ P1Interpolation.py   -> fait le benchmark de la méthode de P1-interpolation
          
          |__ P1InterpolationGaussianProcess.py   -> classe implémentant la méthode
          
          |__ testTool.py   -> fonction pour le calcul de l'erreur L2
          
     ---> SpectralMethod
     
          |__ estimationSpectraleDim1.py   -> classe implémentant l'estimateur de densité spectrale (dimension d'entrée =  1)
          
          |__ estimationSpectraleDim2.py   -> classe implémentant l'estimateur de densité spectrale (dimension d'entrée = 2)
          
          |__ estimationSpectraleDim3.py   -> classe implémentant l'estimateur de densité spectrale (dimension d'entrée = 3)
          
          |__ mySpectralGaussianProcess1D.py   -> classe implémentant la méthode (dimension d'entrée =  1)
          
          |__ mySpectralGaussianProcess2D.py   -> classe implémentant la méthode (dimension d'entrée =  2)
          
          |__ mySpectralGaussianProcess3D.py   -> classe implémentant la méthode (dimension d'entrée =  3)
          
          |__ spectral.py   -> fait le benchmark de la méthode spectrale
          
          |__ spectralModels.py   -> implémentation de plusieurs densités spectrales
          
          |__ testTool.py   -> fonctions utiles pour faire des estimations d'erreur
          
          |__ tool.py   -> fonctions utiles pour construire les pas spatiaux et les pas spectraux évoqués pour cette méthode



 ### rapportdestageTEX  
 
 Rapport de stage sur les méthodes de simulation sous format de fichiers .tex .


     ---> images    -> dossier contenant des images pour le rapport

     ---> Appendice.tex    -> appendice du rapport

     ---> chapitreIntro.tex    -> chapitre d'introduction

     ---> choleskyMethod.tex    -> description de la méthode de Cholesky

     ---> GGGibbsMethod.tex    -> description de la méthode de Galli-Gao-Gibbs

     ---> hmatricesMethod.tex    -> description de la méthode des Hmatrices

     ---> myAbstract.tex    -> abstract du rapport

     ---> P1interpolation.tex    -> description de la méthode de P1-interpolation

     ---> rapportDeStagesubdivided.tex   -> fichier permettant de construire le rapport 

     ---> references.bib    -> références bibliographiques associées au rapport

     ---> spectralMethod.tex    -> description de la méthode spectrale
     

#### Pour produire un pdf du rapport en ligne de commandes:
 - pdflatex rapportDeStagesubdivided.tex
 - bibtex rapportDeStagesubdivided.aux
 - pdflatex rapportDeStagesubdivided.tex


### createTheGraphicsInTheReport

Scripts python permettant de produire les graphiques présents dans le rapport de stage pour chaque méthode de simulation

     ---> choleskyGraphics
     
          |__ choleskyGraphics.py  ->  créer les graphiques pour la méthode de Cholesky  
          
     ---> GGGibbsGraphics
     
          |__ GGGibbsGraphics.py  ->   créer les graphiques pour la méthode de Galli-Gao-Gibbs
          
          
     ---> HmatrixGraphics
     
          |__ HmatrixGraphics.py   ->  créer les graphiques pour la méthode des H-matrices
          
          
     ---> P1interpolationGraphics
     
          |__ P1interpolationGaussianProcess.py   -> classe implémentant la méthode P1 (utile pour créer les graphiques)
          
          |__ P1interpolationGraphics.py   -> créer les graphiques pour la méthode de P1-interpolation
          
          
     ---> spectralGraphics
     
          |__ mySpectralGaussianProcess3D.py   -> classe implémentant la méthode spectrale (dimension d'entrée 3)
          
          |__ spectralGraphics.py  ->  créer les graphiques pour la méthode spectrale
          
          |__ spectralModels.py   -> classe implémentant des densités spectrales

          |__ tool.py   -> fonctions utiles pour le fichier mySpectralGaussianProcess3D.py


#### Pour produire les graphiques en ligne de commandes:

Pour chacun des 5 dossiers présentés ci-dessus, afin d'engendrer les graphiques, il faut se placer dans le dossier
et exécuter le fichier se terminant par -Graphics.py . Après l'exécution, tous les fichiers produits ont le
même nom qu'un fichier apparaissant dans le dossier rapportdestageTEX/images . Certains fichiers peuvent
avoir l'extension .vtk : ce ne sont pas des fichiers images. Pour les fichiers .vtk, il est nécessaire d'utiliser 
un logiciel de visualisation comme paraview permettant de visualiser les maillages ou les réalisations d'un processus.
Dans ce cas, l'utilisateur est contraint de capturer l'image via le logiciel.


### KarhunenLoeveMethod

Description d'une autre méthode de simulation de champs gaussiens à explorer: la méthode de Karhunen-Loeve

          |__ klGaussianProcess.py   -> classe implémentant la méthode de Karhunen-Loeve 
          
          |__ KLSimulationAlgorithm.ipynb   -> notebook expliquant cette méthode
       
          |__ testTool.py  -> fonctions utiles pour d'éventuels tests
          


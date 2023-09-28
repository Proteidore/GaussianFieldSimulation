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
       
          
          


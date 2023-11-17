# Matrix-Matrix-Multiplikation mit Optimierungen 

Aufruf: ./matmul testMatrix1 testMatrix2 OutputFile VAR

Je nach VAR wird eine andere Strategie zur Berechnung der Matrix-Multiplikation gewählt. 
Beim Einlesen der Argumente wird VAR gespeichert und ausgewertet und dann wird die entsprechende Methode zur Berechnung aufgerufen.

-> Es gibt folgende Optionen für VAR: 

### STD - Naiver Algorithmus
Matrix-Multiplikation wird mit Hilfe von 3 FOR-Schleifen durchgeführt


### BLAS - mit BLAS-Library
Benutzt zur Berechnung die Blas-Library-Funktion cblas_gemm


### OPT1 - Transposed
Die 2. Matrix wird erst transponiert, indem die Einträge zeilen- und spaltenweise vertauscht in eine neue Matrix kopiert werden - hierbei verwenden wir noch "Blocking", um beim Kopieren möglichst viele Cache Misses zu vermeiden

Anschließend werden die erste Matrix und die transponierte Matrix naiv multipliziert

### OPT2 - Strassenalgorithmus
Hier unterscheiden wir zunächst einmal, ob die gegebenen Matrizen quadratisch sind und verwenden je nach dem eine andere Methode
Bei beiden Methoden gilt: sind die Matrizen kleiner als eine bestimmte Größe MIN_STRASSEN_SIZE, lohnt sich der Strassen-Algorithmus nicht, daher verwenden wir dann die Naive Implementierung
 
 -> nicht quadratisch: Strassen

    Wir padden die Matrizen auf die nächste quadratische 2^n große Matrix, wobei hier unterschiedlich viele Reihen und Spalten hinzugefügt werden 
    Danach verwenden wir die StrassenQuad-Methode, da der Algorithmus der gleiche ist, sobald die Matrizen eine quadratische, 2^n-Form besitzen

 -> quadratisch: StrassenQuad:

    Die Methode rekursiv aufgerufen um die Matrix in immer kleinere Untermatrizen zu teilen und berechnet mithilfe der 7 Grundadditionen und -multiplikationen die Werte
    Zum Schluss setzen wir die Ergebnismatrix aus den verschiedenen Berechnungen wieder zusammen

### OPT3 - kombiniert Strassen und Transposed
Funktioniert an sich wie der Strassen-Algorithmus, nur dass hier die Transposed-Variante zur restlichen Berechnung verwendet wird, sobald die Matrizen, mit denen wir rechnen eine bestimmte Größe (MIN_STRASSEN_SIZE) unterschreiten


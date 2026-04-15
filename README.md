# TP de Protéomique MADP BIOINFO@LYON
## Contexte Biologique
Vous allez utiliser des outils informatiques qui vous permettront d’analyser des données brutes issues d’une analyse Shotgun Proteomics publiée dans le journal Science sous le titre "Real-time visualization of drug resistance acquisition by horizontal gene transfer reveals a new role for AcrAB-TolC multidrug efflux pump".

Les données associées à cette publication sont publiques et accessibles sur la plateforme [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD011286). Le PDF de la publication est [`data/Nolivos_2019.pdf`](data/Nolivos_2019.pdf).



## Mise en place

### Méthodologie

Vous forkerez le présent "repository" pour vous permettre de sauvegarder votre travail.
Vous le clonerez ensuite dans votre espace de travail.
Vous éditerez ce fichier `README.md` pour répondre aux questions dans les encarts prévus à cet effet et inserer les figures que vous aurez générées. Ce "repository" vous appartenant, vous pouvez créer tous les repertoires et fichiers necessaires à la conduite du TP.

### Installation
La procédure d'installation préférentielle requiert l'utilitaire [uv](https://docs.astral.sh/uv/getting-started/installation/). Si vous ne pouvez pas l'installer, veuillez solliciter votre encadrant.

```sh
git clone https://github.com/<MON_LOGIN_GITHUB>/tp-proteomics.git
cd tp-proteomics
uv sync
uv run jupyter notebook
```

### Test de l'installation

Dans l'interface de jupyter, créez un nouveau fichier notebook (*.ipynb) localisé dans votre repertoire git.
Dans la première cellule copiez le code suivant:

```python
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
```

La première ligne sert à activer le rendu graphique, pour tout le fichier notebook. Pour **dessiner des graphiques**, il vous est conseillé de suivre la méthode illustrée par le code suivant que vous pouvez executer dans une deuxième cellule notebook.

```python
fig, ax = plt.subplots()
x = np.linspace(norm.ppf(0.01),
                norm.ppf(0.99), 100)

ax.plot(x, norm.pdf(x),
       'r-', lw=5, alpha=0.6)

ax.plot(x, np.full(len(x), 0.2),
       'b-', lw=1)

fig.show()
```

* Creation des objets `fig`et `ax`
* Ajout successif de graphiques sur la même figure par l'appel à des methodes de l'objet `ax`
* Affichage de la figure complète via `fig.show()`
* Evaluation de la cellule pour visualisation dans la cellule de résultats.

L'affichage dans la cellule de rendu du notebook devrait confirmer la bonne installation des dépendances.

#### On entend par figure un graphique avec des **axes légendés et un titre**.

La documentation matplotlib est bien faite, mais **attention** il vous est demandé, pour la construction des graphiques, **d'ignorer les méthodes de l'objet `plt`** (frequemment proposées sur le net) et d'utiliser uniquement les méthodes de l'[objet Axes](https://matplotlib.org/api/axes_api.html?highlight=axe#plotting). `plt` est un objet global susceptible d'agir simultanement sur tous les graphiques d'un notebook. A ce stade, son utilisation est donc à éviter.

## Données disponibles

### Mesures experimentales

Un fichier `data/TCL_wt1.tsv` contient les données d'abondances differentielles mesurées sur une souche sauvage d'*Escherichia coli* entre deux conditions: avec Tetracycline et en milieu riche. Le contrôle est le milieu riche.

| Accession | Description | Gene Symbol  |   Corrected Abundance ratio (1.53)    | Log2 Corrected Abundance Ratio | Abundance Ratio Adj. P-Value |   -LOG10 Adj.P-val |
| --- | --- | --- | --- | --- | --- | ---|
| [Accesseur Uniprot](https://www.uniprot.org/help/accession_numbers)  | Texte libre | Texte libre  |  $\frac{ \text{WildType}\_{\text{Tc}} }{ \text{WildType}\_{\text{rich}} }$ |$Log_2(\frac{\text{WildType}\_{\text{Tc}}}{\text{WildType}\_{\text{rich}}})$  | $\mathbb{R}^{+}$ |  $\mathbb{R}^{+}$  |

<!-- | [Accesseur Uniprot](https://www.uniprot.org/help/accession_numbers)  | Texte libre | Texte libre  | <img src="https://render.githubusercontent.com/render/math?math=\frac{\text{WildType}_{\text{Tc}}}{\text{WildType}_{\text{rich}}}"> | <img src="https://render.githubusercontent.com/render/math?math=Log_2(\frac{\text{WildType}_{\text{Tc}}}{\text{WildType}_{\text{rich}}})">  | <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^%2B"> | <img src="https://render.githubusercontent.com/render/math?math=\mathbb{R}^%2B">  | -->


Attention certaines valeurs numériques sont manquantes ou erronées, constatez par vous même en parcourant rapidement le fichier.

### Fiches uniprot

Les fiches de toutes les protéines de *E.coli* sont stockées dans un seul document XML `data/uniprot-proteome_UP000000625.xml`. Nous allons extraire de ce fichier les informations dont nous aurons besoin, à l'aide du module de la librarie standard [XML.etree](https://docs.python.org/3/library/xml.etree.elementtree.html#module-xml.etree.ElementTree). Pour vous facilitez la tâche, les exemples suivants vous sont fournis. Prenez le temps de les executer et de les modifier dans un notebook. Vous vous familliariserez ainsi avec la structure du document XML que vous pouvez egalement inspecter dans un navigateur.

```python
from xml.etree.ElementTree import parse, dump
# Parse the E.Coli proteome XML Document
tree = parse('data/uniprot-proteome_UP000000625.xml')
root = tree.getroot()
ns = '{http://uniprot.org/uniprot}' # MANDATORY PREFIX FOR ANY SEARCH within document
# Store all entries aka proteins in a list of xml nodes
proteins = root.findall(ns + 'entry')
# Display the xml subtree of the first protein 
dump(proteins[0])
```

```python
# Find the xml subtree of a protein with accession "P31224"
for entry in proteins:
    accessions = entry.findall(ns+"accession")
    for acc in accessions:
        if acc.text == "P31224":
            dump(entry)
            break
```

<!-- Cherchez par exemple le sous arbre de la protéine nommée `DACD_ECOLI` -->

### Statistiques de termes GO

Les nombres d'occurences de tous les termes GO trouvés dans le protéome de E.coli sont stockés dans le fichier `data/EColiK12_GOcounts.json`. Ces statistiques ont été dérivées du fichier `data/uniprot-proteome_UP000000625.xml`.

## Objectifs

1. Manipuler les données experimentales tabulées
2. Representer graphiquement les données d'abondance
3. Construire la pvalue des fonctions biologiques (termes GO) portées par les protéines surabondantes
4. Visualiser interactivement les pathways plus enrichis.

### Description statistique des Fold Change
La lecture des données au format tabulé est l'occasion de se familliariser avec la [librairie pandas](https://pandas.pydata.org).

##### Lecture de données

###### source:`data/TCL_wt1.tsv`

La fonction `read_csv` accepte différents [arguments](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html) de format de données très utiles.

```python
df = pandas.read_csv()
```

Quel est le type de l'objet `df`?
```
pandas.DataFrame
```

##### Descriptions d'une table de données
Que permettent les méthodes suivantes?
###### df.shape
```
(1750, 7)
```
###### df.head()
```
5 premières lignes du tableau.
```
###### df.tail()
```
5 dernières lignes du tableau.
```
###### df.columns
```
Nom des colonnes.
```
###### df.dtypes
```
Types des colonnes.
```
###### df.info
```
Résumé du DataFrame.
```
###### df.describe()
```
Affiche des statistiques descriptives de la colonne numériques.
```
###### df.dropna()
```
Enlève les valeurs manquantes.
```

##### Accès aux éléments d'une table de données

```python
values = df[['Description', 'Gene Symbol']]
```

Quel est le type de `values` ?

Verifiez si certaines méthodes de `DataFrame` lui sont applicables.
Ce type supporte l'accès par indice et les slice `[a:b]`

##### Accès indicé

On peut accéder aux valeurs du DataFrame via des indices ou plages d'indice. La structure se comporte alors comme une matrice. La cellule en haut et à gauche est de coordonnées (0,0).
Il y a différentes manières de le faire, l'utilisation de `.iloc[slice_ligne,slice_colonne]` constitue une des solutions les plus simples. N'oublions pas que shape permet d'obtenir les dimensions (lignes et colonnes) du DataFrame.
###### Acceder aux cinq premières lignes de toutes les colonnes
```python
values.head()

```

###### Acceder à toutes les lignes de la dernière colonne
```python
values.iloc[:,-1]

```

###### Acceder aux cinq premières lignes des colonnes 0, 2 et 3
```python
df.head().iloc[:,[0,2,3]]

```

##### Conversion de type

Le type des valeurs d'une colonne peut être spécifiée:

* à la lecture

```python
pandas.read_csv('data/TCL_wt1.tsv', sep="\t",  dtype = {'Accession': str, 'Description': str, 'Gene Symbol': str, 
                                                 'Corrected Abundance ratio (1.53)': float,  'Log2 Corrected Abundance Ratio': float, 
                                                 'Abundance Ratio Adj. P-Value: (127. T3 Tc WT) / (126. T0 WT)': float, '-LOG10 Adj.P-val': float})
```

* modifiée à la volée

```python
df = df.astype({'Log2 Corrected Abundance Ratio': float, '-LOG10 Adj.P-val': float } )
```

On peut specifier la chaine de caractères correspondant au valeurs "incorrectes" produites par excel comme-ci

```python
pandas.read_csv('data/TCL_wt1.tsv', sep="\t", na_values="#VALEUR!")
```
##### Selection avec contraintes
La méthode `loc` permet de selectionner toutes les lignes/colonnes respectant certaines contraintes

* Contraintes de valeurs continues

```python
df.loc[(df['-LOG10 Adj.P-val'] > 0 )  & (df['Log2 Corrected Abundance Ratio'] > 0.0 ) ]
```

* Contraintes de valeurs discrètes

```python
df.loc[ df['Gene Symbol'].isin(['fadR', 'arcA'] ) ]
```


#### Appliquons ces outils à l'analyse de données protéomique

##### 1. Chargez le contenu du fichier `data/TCL_wt1.tsv` dans un notebook en eliminant les lignes porteuses de valeurs numériques aberrantes

##### 2. Representez par un histogramme les valeurs de `Log2 Corrected Abundance Ratio`

<!-- ##### 3. A partir de cette échantillon de ratio d'abondance,  estimez la moyenne <img src="https://render.githubusercontent.com/render/math?math=\mu"> et l'ecart-type <img src="https://render.githubusercontent.com/render/math?math=\sigma"> d'une loi normale. -->

##### 3. A partir de cette échantillon de ratio d'abondance,  estimez la moyenne $\mu$ et l'ecart-type $\sigma$ d'une loi normale.
```
mu = -0.6467130248461945
sigma = 0.46723442417098815
```

##### 4. Superposez la densité de probabilité de cette loi sur l'histogramme. Attention, la densité de probabilité devra être mis à l'echelle de l'histogramme (cf ci-dessous)


```python
# _ est le vecteur des valeurs d'abondance
fig, ax = plt.subplots()
hist = ax.hist(_, bins=100) # draw histogram
x = np.linspace(min(_), max(_), 100) # generate PDF domain points
dx = hist[1][1] - hist[1][0] # Get single value bar height
scale = len(_)*dx # scale accordingly
ax.plot(x, norm.pdf(x, mu, sigma)*scale) # compute theoritical PDF and draw it
```

![Histogramme des Log2 FC](figures/histogram_log2FC.png)

##### 5. Quelles remarques peut-on faire à l'observation de l'histogramme et de la loi théorique?

```
La distribution des valeurs ne suit pas une loi Normale.
```

#### Construction d'un volcano plot

##### A l'aide de la méthode [scatter](https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.scatter.html) representer $\text{Log}_{10}({\text{p-value}}) = f(\text{Log}_2(\text{abundance ratio}))$

##### Matérialisez le quadrant des protéines surabondantes, par deux droites ou un rectangle
Sont condidérées comme surabondantes les proteines remplissant ces deux critères:

* $\text{Log}_2(\text{abundance ratio})\gt\mu%2B\sigma$
* $\text{p-value}<0.001$

![Volcano plot](figures/volcano_plot.png)

### Analyse Fonctionelle de pathway

Nous allons implementer une approche ORA (Over Representation Analysis) naive.

##### 1. Retrouvez les entrées du fichier TSV des protéines surabondantes

Quelles sont leurs identifiants UNIPROT ?
``` 
['P23721',
 'P77804',
 'P0A6K6',
 'P0A799',
 'P0A7G6',
 'P0A6F3',
 'P25745',
 'P0A6M8',
 'P0A6L0',
 'P0A8V6',
 'P0A9Q1',
 'P02358',
 'P0ACF8',
 'P62399',
 'P0A905',
 'P76506',
 'P13036',
 'P10384',
 'P06971',
 'P0A910'
 'P06996',
 'P76344',
 'P02931']
```

#### 2. Listez les termes GO portés par ces protéines surabondates

Les `entry` du fichier `data/uniprot-proteome_UP000000625.xml` présentent des balises de ce type:

```xml
<dbReference type="GO" id="GO:0005737">
<property type="term" value="C:cytoplasm"/>
<property type="evidence" value="ECO:0000501"/>
<property type="project" value="UniProtKB-SubCell"/>
</dbReference>
```

Pour vous aider à obtenir la liste des termes GO des protéines surabondantes, voici la fonction `getAccessionGOTerms(xml_file, uniprotID)`
qui extrait du fichier de protéome uniprot xml fourni en 1er argument les termes GO portés par la protéine au code uniprot fourni en deuxième argument.

```python
from xml.etree.ElementTree import parse, dump

def getAccessionGOTerms(xmlFile, accession):
    tree = parse(xmlFile)
    root = tree.getroot()
    ns = '{http://uniprot.org/uniprot}'
    
    match_go_terms = []
    proteins = root.findall(ns + 'entry')
    for entry in proteins:
        accessions = entry.findall(ns+"accession")
        current_accessions = [ acc.text for acc in accessions ]
        if not accession in current_accessions:
            continue
        goTerms = entry.findall('{http://uniprot.org/uniprot}dbReference[@type="GO"]')
        #goTerms = xmlEntry.findall(ns +'dbReference[@type="GO"]')
        for goT in goTerms:
            gID   = goT.attrib['id']
            gName = goT.find(ns +'property[@type="term"]').attrib['value']
            match_go_terms.append((gID, gName))
        break
    return match_go_terms
getAccessionGOTerms("./data/uniprot-proteome_UP000000625.xml", "P0A8V6")
```
A l'aide de cette fonction, il devrait être possible de construire le dictonnaire des termes GO de toutes les protéines surabondantes précedemment identifiées.
Ce dictionnaire pourrait être de la forme suivante:
```python
{'GO:0005829': {'ID'        : 'GO:0005829',
                'name'      : 'C:cytosol',
                'carried_by': ['P0A8V6', 'P0A9Q1', 'P02358', 'P0ACF8', 'P62399', 'P76344']
                },
 'GO:0003677': {'ID'        : 'GO:0003677',
                'name'      : 'F:DNA binding',
                'carried_by': ['P0A8V6', 'P0ACF8']
                }
  }
```
Vous implémenterez la construction de ce dictionnaire et ainsi stockerez, pour la suite de l'analyse, les représentations des termes GO parmi les protéines surabondantes.

#### 3. Obtention des paramètres du modèle

Nous évaluerons la significativité de la présence de tous les termes GO portés par les protéines surabondantes à l'aide d'un [modèle hypergéometrique](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html).

Si k protéines surabondantes porte un terme GO, la pvalue de ce terme sera équivalente à  $P(X\ge k | X \sim H(k,K,n,N) )$

Completer le tableau ci-dessous avec les quantités vous semblant adéquates pour modeliser la pvalue de **chaque pathway [termes GO]**

| Symboles | Paramètres                   | Quantités Biologiques                                |
| -------- | ---------------------------- | ---------------------------------------------------- |
| k        | nombre de succès observés    | nombre de protéines surabondante portant le terme GO |
| K        | nombre de succès possibles   | nombre total de protéines portant le terme GO        |
| n        | nombre d'observations        | nombre total de protéines surabondantes              |
| N        | nombre d'elements observables| nombre total de protéines mesurées                   |

#### 4. Calcul de l'enrichissement en fonctions biologiques

A l'aide du contenu de `data/EColiK12_GOcounts.json` parametrez la loi hypergeometrique et calculez la pvalue
de chaque terme GO portés par les protéines surabondantes. Vous reporterez ces données dans le tableau ci-dessous

| identifiant GO | définition                                              | occurence | pvalue   |
|----------------|---------------------------------------------------------|-----------|----------|
| GO:0005829     | C:cytosol                                               | 13        | 9.32e-17 |
| GO:0009279     | C:cell outer membrane                                   | 8         | 1.53e-15 |
| GO:0006974     | P:cellular response to DNA damage stimulus              | 6         | 1.69e-09 |
| GO:0034220     | P:ion transmembrane transport                           | 3         | 4.95e-08 |
| GO:0042802     | F:identical protein binding                             | 6         | 1.27e-07 |
| GO:0046930     | C:pore complex                                          | 3         | 2.18e-07 |
| GO:0015288     | F:porin activity                                        | 3         | 2.56e-07 |
| GO:0009264     | P:deoxyribonucleotide catabolic process                 | 2         | 4.01e-07 |
| GO:0038023     | F:signaling receptor activity                           | 2         | 6.01e-06 |
| GO:0015344     | F:siderophore uptake transmembrane transporter activity | 2         | 1.12e-05 |
|                |                                                         |           |          |

Quelle interpretation biologique faites-vous de cet enrichissement en termes GO ?
Aprés analyse des termes GO, on observe une surreprésentation des termes associés à la membrane externe bactérienne, complexes de pores et au transport transmembranaire. Ces fonctions jouent un rôle dans la perméabilité membraneire et aux systèmes de transport moléculaire.
Cela peut indiquer que des processus s'activent pour expulser les antibiotique hors de la cellule et réduire leur concentretion intracellulaire. Ce qui est cohérent dans le cadre de notre étude d'efflux AcrAB-TolC.

### Analyse des interactions répertoriées dans STRING

#### Chargement des données dans STRING

Dans le menu de gauche, sélectionner 'Multiple proteins'.

Copier-coller les identifiants Uniprot des 
protéines sur-exprimées et spécifier l'organisme 'Escherichia coli K-12".

Valider le mapping produit par STRING en clickant sur 'Continue'.


#### Visualisation du réseau dans STRING


Combien d'interactions contient ce réseau ?

```

38


```


Faire varier les paramètres de visualisation du réseau dans 'Settings' pour afficher le réseau fonctionnel
ou physique avec différents indices de confiance.


Combien d'interactions sont supportées par chaque source ('Textmining', 'Experiments', 'Databases','Co-expression',
'Neighborhood', 'Gene Fusion', 'Co-occurence') ?

Hint: l'onglet Analysis, donne accès aux nombre des interactions du réseau.
```
'Textmining':29
'Experiments':8
'Databases': 4
'Co-expression': 9
'Neighborhood':1
'Gene Fusion': 0
'Co-occurence': 7

```

#### Analyse du réseau des protéines sur-exprimées dans le contexte du réseau global.

Consulter la rubrique 'Network Stats' dans l'onglet Analysis. 

Que peut-on en conclure sur les interactions de ce petit ensemble de protéines ?
```


on a significativement plus d'interaction qu'attendu.

```

Afin de replacer ces protéines dans le contexte du réseau d'interaction global de E. coli, 
ajouter les interacteurs de la première et de la deuxième couche.

Que pouvez-vous en déduire sur les mécanismes activés par la présente de tétracycline ?
```




```
#### Analyse de sur-représentation des termes GO

Consulter l'analyse de sur-représention des termes GO présents dans l'onglet 'Analysis'.
Est-ce cohérent avec votre analyse précédente ?

```




```



### Construction de réseaux d'interactions à partir de données MITAB.

Le format (MITAB) stocke des paires de protéines en interaction. Dans ce format, chaque colonne porte une information spécifique.
Une description du format est disponible [ici](https://psicquic.github.io/MITAB27Format.html).
Les données d'interactions impliquant les protéines surreprésentées de l'expérience ont été obtenues depuis la base de données [Intact](https://www.ebi.ac.uk/intact/home).
Ces données sont mises à votre disposition dans le fichier `data/proteins.mitab`.

Vous extrairez du fichier les paires d'identifiants uniprot des protéines en interaction.

Ces paires de protéines constituent un réseau d'interaction protéine-protéine que vous allez dessiner à l'aide de la libraire [networkx](https://networkx.org/documentation/stable/reference).
Le code suivant vous est fourni à titre d'exemple.

```python
%matplotlib inline
import matplotlib.pyplot as plt
import networkx as nx
G = nx.Graph()

fig, ax = plt.subplots(figsize=(8, 8))

G.add_edge('a', 'b')
G.add_edge('e', 'b')
G.add_edge('e', 'a')
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_color=['blue','blue','red'] , node_size=2000)
```

![exemple reseau](data/reseau_exemple.jpg)

Les positions des noeuds sont paramétrables au travers de l'objet [layout](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html). Une fois une première représentation du réseau obtenue, affinez celle-ci afin de:

* colorier dans une couleur spécifique uniquement les protéines surabondantes dans l'expérience.
* Écrire les identifiants uniprot dans les noeuds du réseau.
* colorier les protéines appartenant à des classes GO communes.
* Faire du diamètre des noeuds une fonction du nombre de partenaires protéiques.
* N'afficher que les noeuds des protéines mesurées dans l'experience
* Utiliser une échelle de couleur continue fonction de l'abondance pour colorier les noeuds

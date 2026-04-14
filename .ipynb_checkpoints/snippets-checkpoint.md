```python
import re
def mitab_parser(mitab:str, uniprot_ac:list[str])->list[ list[str] ]:
    """ Extract from provided file (path), lines featurin at least one interactor 
        present in uniprot_ac provided list
        Returns: A list of mitab records split over tabulation.
    """
    hits = []
    with open(mitab, 'r') as fp:    
        for line in fp:
            line = line.split()
            for proteinfield in line[0:2]: # Test column One and Two
                _ = re.match("uniprotkb:(\\S+)", proteinfield) # Extract uniprot identifier
                x = _[1] if _ else None
                if x in uniprot_ac: # Check uniprot identifier is in provided ac list
                    hits.append(line)
                    break
                    
    print(f"Filtered {len(hits)} interactions")
    return hits

res = mitab_parser("/Users/glaunay/t.txt", ['P0A6F5', 'P02931'])
```

```python
import pandas
df = pandas.read_csv('data/TCL_wt1.tsv', sep="\t",  na_values="#VALEUR!")
df = df.dropna()
df_abd = df.loc[(df['-LOG10 Adj.P-val'] > 3 )  & 
    (df['Log2 Corrected Abundance Ratio'] > -0.6 ) ]
abd_acc = df_abd['Accession'].tolist()
```

```python
GO_dic = {}
for acc in abd_acc:    
    go_terms = getAccessionGOTerms(
        "./data/uniprot-proteome_UP000000625.xml", acc)
    for go_id, go_name in go_terms:
        if not go_id in GO_dic:
            GO_dic[go_id] = { 
                'ID' : go_id, 'name' : go_name, 'carried_by' : []
            }
        GO_dic[go_id]['carried_by'].append(acc)
```

```python
# Background statistics
import json
from scipy.stats import hypergeom
background = {}
with open("data/EColiK12_GOcounts.json", "r") as fp:
    background = json.load(fp)["go_terms"]
"""eg:
 {'GO:0000006': {'count': 1,
   'name': 'F:high-affinity zinc transmembrane transporter activity'}

"""
def compute_pvalue(GO_term:dict, background, ttl_abdnt :int, N = 1800):
    """ Compute pvalue of enrichissment of passed GO_term
    Where GO_term is of the shape :
        { 'ID' : 'GO:0000006', 'name' : "blabla", 'carried_by' : ["P0000", "P11111"] }
    """
    
    GO_id = GO_term['ID']
    GO_name = GO_term['name']
    
    n = ttl_abdnt # Total nb of draw
    k_obs = len(GO_term["carried_by"]) # observed nb of success
    K = background[GO_id]['count']
    rv = hypergeom(N, K, n)
    p_value = 0 # pvalue = P(X>=k_obs)
    for k in range(k_obs, ttl_abdnt + 1):
        p_value += rv.pmf(k)
    return p_value, GO_id, GO_name
    
```

```python
compute_pvalue(GO_dic['GO:0000062'], background, 23, 1800)
# Applied to all GO go_terms
all_scores = []
for GO_ID in GO_dic:
    all_scores.append( 
        compute_pvalue(GO_dic[GO_ID], background, 23, 1800)
    )
sorted(all_scores, key = lambda t:t[0])
```

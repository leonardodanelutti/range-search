#import "template.typ": *
#import "@preview/lovelace:0.3.0": *

#show: ams-article.with(
  paper-size: "a4",
  title: [ORTHOGONAL RANGE SEARCH],
  authors: (
    (
      name: "Leonardo Danelutti",
      organization: [Università di Udine],
      email: "danelutti.leonardo@spes.uniud.it",
    ),
  ),
  bibliography: bibliography("refs.bib"),
)

#set text(lang: "it")

= Introduzione
La ricerca ortogonale (o di intervalli rettangolari) è il problema di trovare tutti i punti di un insieme $P$ di dimensione fissa la cui posizione cade entro un intervallo di coordinate definito da un rettangolo (in dimensione 2) o un iperrettangolo (in dimensioni superiori). Un tipico esempio è una query su un database: ogni record diventa un punto con coordinate (data di nascita, salario, ecc.) e si cerca, ad esempio, "tutti i dipendenti nati fra il 1950 e il 1955 che guadagnano fra $3000 e $4000 al mese". In generale, dati $n$ punti in $R^d$ (uno per ciascun record con $d$ attributi), una query ortogonale restituisce tutti i punti che ricadono in un intervallo di coordinate del tipo $[x_1:x'_1] times [x_2:x'_2] times dots times[x_d:x'_d]$.
#footnote[
  Questo elaborato segue quanto riportato nel capitolo 5 del libro di de Berg et al. @computational_geom, anche le figure sono tratte da lì.
]

#figure(
  image("imgs/RangeQuery.png", width: 60%),
  caption: [Esempio di query ortogonale in un database con due attributi: data di nascita e salario.]
)



= Similarità con gli Alberi Binari di Ricerca
Per comprendere le strutture dati per la ricerca ortogonale, è utile partire dal caso più semplice: la ricerca di intervallo unidimensionale. Dato un insieme di valori numerici e un intervallo $[a, b]$, vogliamo trovare tutti i valori che cadono nell'intervallo.

Un albero binario di ricerca bilanciato risolve questo problema in modo elegante: per trovare tutti i valori in $[a, b]$, si individua il nodo $v$ che rappresenta il "lowest common ancestor" dei valori $a$ e $b$ ($v_"split"$ nella @bst-range-search). Da questo punto, si esplorano tutti i sottoalberi che intersecano l'intervallo: nel sottoalbero sinistro si cercano i valori $>= a$, nel sottoalbero destro quelli $<= b$. Il tempo di query risulta $O(log n + k)$ dove $k$ è il numero di elementi trovati, dato che la ricerca dei nodi $mu$ e $mu'$ richiede $O(log n)$ e la visita dei sottoalberi che contengono i valori nell'intervallo richiede tempo lineare nel numero di nodi riportati.

#figure(
  image("imgs/BSTRangeSearch.png", width: 50%),
  caption: [Esempio di ricerca di intervallo unidimensionale in un albero binario di ricerca.]
) <bst-range-search>

= KD-Tree
Il KD-Tree rappresenta la generalizzazione più diretta e intuitiva degli alberi binari di ricerca al caso multi-dimensionale. L'idea è semplice: se nel caso unidimensionale partizioniamo lo spazio rispetto ad una singola dimensione, nel caso multi-dimensionale possiamo partizionare lo spazio prima rispetto ad una posizione, poi ad un'altra, e così via, ciclicamente.

Un kd-tree (k-dimensional tree) è una struttura dati binaria che suddivide ricorsivamente lo spazio alternando dimensioni di split. In 2D, ad ogni nodo interno si seleziona una linea verticale (split su $x$) o orizzontale (split su $y$) in modo alternato ai livelli. Questo approccio mantiene la semplicità concettuale degli alberi binari: ogni nodo divide ancora lo spazio in due parti, ma la dimensione su cui avviene la divisione cambia ciclicamente.

#figure(
  image("imgs/KDTreeConstruction.png", width: 100%),
  caption: [Esempio di KD-Tree in 2D con partizioni alternate su $x$ e $y$.]
) <kd-tree-construction>

== Creazione di un KD-Tree

Ogni nodo interno memorizza il valore di split e la dimensione associata, mentre le foglie contengono singoli punti di $P$. Ecco uno schema dell'algoritmo di costruzione:

#figure(
  kind: "algorithm",
  supplement: [Algoritmo],

  pseudocode-list(booktabs: true, stroke:none, numbered-title: [Costruzione di un KD-Tree])[
  + *function* BUILD-KD-TREE(P, depth):
    +  *if* |P| = 1:
      + crea un nodo foglia che contiene l'unico punto di P e restituiscilo.
    + *else*:
      + dim = depth mod k
      + _Trova il valore mediano_ m _di_ P _nella coordinata _dim
      + _Suddividi_ P _in_ P_left = { p | p.coord[dim] < m } _e_ P_right = { p | p.coord[dim] >= m }.
      + costr_left  = BUILD-KD-TREE(P_left, depth+1)
      + costr_right = BUILD-KD-TREE(P_right, depth+1)
      + _crea un nodo interno_ v _con_ split = m, dim = dim, _figli_ v.left=costr_left _e_ v.right=costr_right.
      + *return* v.
  ]
)

La complessità per trovare il valore mediano e dividere i punti è $O(n)$, l'equazione di ricorrenza è quindi $T(n) = 2T(n/2) + O(n)$. La costruzione del kd-tree richiede quindi $O(n log n)$ tempo. Lo spazio utilizzato è $O(n)$ in quanto ogni punto risiede in una singola foglia ed il numero di nodi interni è nell'ordine del numero di foglie (in quanto bilanciato).

== Query in un KD-Tree

Per eseguire una query ortogonale con un kd-tree, si procede in modo simile alla ricerca di intervallo unidimensionale, ma tenendo conto delle diverse dimensioni. L'algoritmo visita i nodi dell'albero e decide, in base alla regione di spazio associata al nodo e all'intervallo di query, se esplorare i figli sinistro e/o destro.

#figure(
    kind: "algorithm",
    supplement: [Algoritmo],

    pseudocode-list(booktabs: true, stroke:none, numbered-title: [Query in un KD-Tree])[
    + *function* KD-RANGE-SEARCH(v, range, results):
      +  *if* v _è una foglia_:
        + *if* v.point _è in_ range:
          + aggiungi v.point a results
      + *else*:
        + dim = v.dim
        + *if* range.min[dim] < v.split:
          + KD-RANGE-SEARCH(v.left, range, results)
        + *if* range.max[dim] >= v.split:
          + KD-RANGE-SEARCH(v.right, range, results)
  ]
)

#figure(
  image("imgs/KDTreeQuery.png", width: 80%),
  caption: [Una query in un KDTree.]
) <kd-tree-query>

#theorem(
  [La complessità della query in un kd-tree è $O(n^(1-1/d) + k)$, dove $d$ è la dimensione dello spazio e $k$ il numero di punti riportati.]
)

Per ogni nodo $v$ interno dell'albero, chiamiamo $"region"(v)$ il rettangolo, che può essere illimitato in una o più direzioni, delimitato dalle line di split degli antenati di $v$. Se $r$ è la radice dell'albero, allora $"region"(r) = R^d$. Ad esempio nella @kd-tree-construction, $"region"(l_4)$ è la regione del piano sotto a la linea $l_2$ (in quanto $l_4$ si trova nel sottoalbero sinistro di $l_2$) e a sinistra di $l_1$ (in quanto $l_4$ si trova nel sottoalbero sinistro di $l_1$).

Durante la ricerca si possono incontrare tre tipi di nodi $v$:
1. nodi in cui $"region"(v)$ non interseca la regione di query: in questo caso non si esplorano i sottoalberi radicati in $v$. (nodi bianchi nella @kd-tree-query)
2. nodi in cui $"region"(v)$ è completamente contenuta nella regione di query. La complessità per riportare tutti i punti all'interno del sottoalbero radicato in $v$ è lineare nel numero dei nodi dell'albero. La complessità di tutte le chiamate su questi nodi è quindi $O(k)$. (nodi grigio scuri nella @kd-tree-query)
3. nodi in cui $"region"(v)$ è parzialmente contenuta nella regione di query. Ogni nodo di questo tipo corrisponde quindi ad una regione intersecata dal bordo della regione di query. (nodi grigio chiari nella @kd-tree-query)

Per trovare un limite superiore al numero di nodi del terzo tipo possiamo contare il numero di regioni intersecate da una retta parallela ad un asse (iperpiani nel caso di dimensione maggiori di due), questo corrisponde al caso pessimo dove il bordo della regione di query interseca il maggior numero di regioni. Chiamiamo questo numero $Q(n)$.
Per calcolare $Q(n)$ costruiamo l'equazione di ricorsione sui nodi che dividono il piano (iperspazio) nella stessa dimensione.
Partendo dalla radice suddividiamo il piano (iperspazio) nelle regioni associate ai nodi fino a profondità $d$, così da poter arrivare ad un nodo con la stessa dimensione, si avranno quindi $2^d$ regioni. Dato che le linee di split sono perpendicolari tra loro $l$ attraversa esattamente $2^(d-1)$ di queste regioni. Si ha quindi che:

$ Q(n) := cases(
  O(1) &"  se" n = 1,
  2^(d-1) Q(n/2^d) + 2^(d-1) &"  se" n > 1
) $

$ 


Q(1) &= O(1) \
Q(n) &= 2^(d-1) Q(n/2^d) + 2^(d-1)
$

Che ha soluzione $Q(n) = O(n^(1-1/d))$.

Il numero totale di nodi visitati è quindi $O(n^(1-1/d) + k)$. Si noti che questo è il caso peggiore, in pratica il numero di regioni che intersecano il bordo della regione di query è spesso molto più basso. Ad esempio se cerchiamo un singolo punto (cioè nel range $[x_1:x_1] times [y_1:y_1]$) il tempo di query diventa $O(log n)$.

#figure(
  kind: image,
  table(
    stroke: 0pt,
    columns: 2,
    image("imgs/KDQueryCompl1.png"),
    image("imgs/KDQueryCompl2.png"),
  ),
  caption: [Esempio con d=2. Il piano è diviso in 4 regioni, la retta verticale (in rosso) interseca le due regioni alla sinistra della radice]
)

= Range Tree

Un range tree è un albero multi-livello che, rispetto ai kd-tree, consente query più veloci a costo di più spazio. Anche in questo caso la struttura è una generalizzazione degli alberi binari di ricerca, che in questo caso contiene informazione ridondante per velocizzare le query. Per il caso in due dimensioni viene costruito l'albero binario di ricerca sulla prima dimensione, ad ogni nodo interno $v$ si associa una struttura di supporto: un albero di ricerca bilanciato sui valori di $y$ dei punti presenti nelle fogli dell'albero radicato in $v$. Ogni nodo interno contiene: il valore di split $x_med$, un puntatore al sottoalbero sinistro e destro, ed un puntatore alla radice di un associato (un albero di ricerca su $y$). Per più dimensioni si procede in modo analogo, costruendo ad ogni nodo un albero di ricerca sulla dimensione successiva.


#table(
  columns: 2,
  stroke: 0pt,
  figure(
    image("imgs/RangeTreeConstruction.png", width: 60%),
    caption: [Ogni nodo dell'albero memorizza un puntatore alla sua struttura associata.]
  ),
  figure(
    image("imgs/RangeTreeAssoc.png", width: 70%),
    caption: [Un Range Tree con dimensione maggiore di 2. Le strutture associate sono annidate.]
  )
)

== Costruzione del Range Tree

#figure(
  kind: "algorithm",
  supplement: [Algoritmo],

  pseudocode-list(booktabs: true, stroke:none, numbered-title: [Costruzione di un Range Tree])[
  + *function* BUILD-RANGE-TREE(P, d):
    +  *if* |P| = 1:
      + _crea una foglia che contiene l'unico punto di_ P _e restituiscila_.
    + *else*:
      + _trova il valore mediano_ m _di_ P _nella coordinata_ d;
      + _suddividi_ P _in_ P_left = { p | p.coord[d] < m } _e_ P_right = { p | p.coord[d] >= m };
      + costr_left  = BUILD-RANGE-TREE(P_left, d);
      + costr_right = BUILD-RANGE-TREE(P_right, d);
      + assoc = BUILD-RANGE-TREE(P, d+1) _se_ d+1 <= D _altrimenti_ BUILD-SEARCH-TREE(P);
      + _crea un nodo interno_ v _con_ split = m, left=costr_left, right=costr_right, assoc=assoc;
      + *return* v.
  ]
)


Il sorgente mostra questo comportamento: ad ogni chiamata ricorsiva si costruisce un albero di supporto $T_"assoc"$ sui punti correnti.

L'equazione di ricorrenza il tempo di costruzione è:

$
T(n, d) := cases(
  O(1) &"  se" n = 1,
  O(n log n) &"  se" d = 1,
  2T(n/2, d) + T(n, d-1) + O(n) &"  se" n > 1 and d > 1
)
$

Dove O(n log n) è il costo di costruzione dell'albero di ricerca bilanciato nel caso base (d=1), $2T(n/2, d)$ è il costo per costruire i sottoalberi sinistro e destro, $T(n, d-1)$ è il costo per costruire l'albero associato, e $O(n)$ è il costo per trovare il mediano e partizionare i punti.

La soluzione all'equazione di ricorsione è $T(n, d) = O(n log^(d-1) n)$.

In modo analogo lo spazio occupato è $O(n log^(d-1) n)$ (cambia solo il caso base della ricorrenza).

== Query in un Range Tree

L'idea della query in un range tree è simile a quella del caso unidimensionale, quando si individua un sottoalbero che è contenuto completamente nella regione di query nella prima dimensione, invece di riportare tutti i punti del sottoalbero, si esegue una query sull'albero associato della prossima dimensione.

#figure(
  kind: "algorithm",
  supplement: [Algoritmo],

  pseudocode-list(booktabs: true, stroke:none, numbered-title: [Query in un Range Tree])[
  + *function* RANGE-TREE-RANGE-SEARCH(v, range, results):
    +  *if* v _è una foglia_:
      + *if* v.point _è in_ range:
        + aggiungi v.point a results
    + *else*:
      + _trova il nodo di split_ v_split

      + \// Segui il percorso verso x_n
      + v = v_split.left
      + *while* v _non è una foglia_:
        + dim = v.dim
        + *if* range.min[dim] < v.split:
          + \// Il sottoalbero destro è completamente contenuto nella regione di query
          + *if* dim == d - 2:
            + BST-RANGE-SEARCH(v.right.assoc, range, results)
          + *else*:
            + RANGE-TREE-RANGE-SEARCH(v.right.assoc, range, results)
          + v = v.left
        + *else*:
          + v = v.right

      + \// v è una foglia
      + *if* v.point _è in_ range:
        + aggiungi v.point a results
      + \
      + \// Segui il percorso verso x'\_n
      + \// ...
  ]
)

La visita dalla radice ai nodi foglia nel range tree d-dimensionale richiede $O(log n)$ tempo. Durante la ricerca possiamo visitare fino a $O(log n)$ nodi dove viene chiamata la ricerca sull'albero associato. La ricorsione che descrive il tempo per una query, senza contare i punti riportati, è:

$
T(n, d) := cases(
  O(1) &"  se" n = 1,
  O(log n) &"  se" d = 1 "(Range query in un BST)",
  O(log n) + O(log n) T(n, d-1) &"  se" n > 1 and d > 1
)
$

\


La ricorrenza ha soluzione $T(n, d) = O(log^d n)$. Aggiungendo il costo di riportare i $k$ punti trovati si ottiene la complessità finale:

$ O(log^d n + k) $


= Layered Range Tree

Utilizzando una tecnica chiamata fractional cascading è possibile ridurre il tempo di query del range tree a $O(log^(d-1) n + k)$. Questa tecnica richiede una modifica della struttura dati associata nell'ultima dimensione. Invece di utilizzare un albero di ricerca bilanciato per memorizzare il sottoinsieme canonico di un nodo, si utilizza un array ordinato dei punti in quella dimensione. Inoltre per ogni elemento dell'array si tiene traccia di due puntatori (left e right) che indicano la posizione del primo elemento nell'array del figlio sinistro e destro che è maggiore o uguale all'elemento corrente. In questo modo, quando si esegue una ricerca binaria nell'array associato di un nodo, si possono seguire i puntatori per trovare rapidamente la posizione corrispondente negli array dei figli, evitando ricerche binarie ripetute.

#figure(
    image("imgs/RangeTree2D.png", width: 90%),
)

#figure(
    image("imgs/CascadeArrays.png", width: 90%),
    caption: [In alto un range tree in 2D, in basso gli array associati per ogni nodo dell'albero.]
)

== Costruzione del Layered Range Tree

#figure(
    kind: "algorithm",
    supplement: [Algoritmo],

    pseudocode-list(booktabs: true, stroke:none, numbered-title: [Costruzione di un Layered Range Tree])[
    + *function* BUILD-LAYERED-RANGE-TREE(P, d):
      +  *if* |P| = 1:
        + _crea una foglia che contiene l'unico punto di_ P _e restituiscila_.
      + *else*:
        + _trova il valore mediano_ m _di_ P _nella coordinata_ d;
        + _suddividi_ P _in_ P_left = { p | p.coord[d] < m } _e_ P_right = { p | p.coord[d] >= m };
        + costr_left  = BUILD-LAYERED-RANGE-TREE(P_left, d);
        + costr_right = BUILD-LAYERED-RANGE-TREE(P_right, d);
        + assoc = BUILD-LAYERED-RANGE-TREE(P, d+1) _se_ d+1 <= D _altrimenti_ BUILD-CASCADE-ARRAY(costr_left.assoc, costr_right.assoc);
        + _crea un nodo interno_ v _con_ split = m, left=costr_left, right=costr_right, assoc=assoc;
        + *return* v.
        + \

    + *function* BUILD-CASCADE-ARRAY(left_array, right_array):
      +  i, j = 0
      +  cascade_array = []
      +  *while* i < |left_array| *and* j < |right_array|:
        + *if* left_array[i].y < right_array[j].y:
          + elem = Element(point=left_array[i], left=i, right=j)
          + cascade_array.append(elem)
          + i += 1
        + *else*:
          + elem = Element(point=right_array[j], left=i, right=j)
          + cascade_array.append(elem)
          + j += 1
      +  *while* i < |left_array|:
        + elem = Element(point=left_array[i], left=i, right=j)
        + cascade_array.append(elem)
        + i += 1
      +  *while* j < |right_array|:
        + elem = Element(point=right_array[j], left=i, right=j)
          + cascade_array.append(elem)
          + j += 1
      + *return* cascade_array
    ]
)


Si noti come l'algoritmo di costruzione del layered range tree è identico a quello del range tree, con la differenza che nell'ultima dimensione viene costruita la nuova struttura dati invece di un albero di ricerca bilanciato. La costruzione dell'array associato ad un nodo è il merge delle due strutture associate dei figli, tenendo inoltre traccia della posizione dei puntatori al momento della creazione di un elemento (Il codice riportato non tiene conto di punti con la stessa coordinata, situazione facilmente risolvibile).

La costruzione dell'array associato richiede $O(n)$ tempo, l'equazione di ricorrenza rimane quindi invariata, se non per il caso base:

$
T(n, d) := cases(
  O(1) &"  se" n = 1,
  O(n) &"  se" d = 1,
  2T(n/2, d) + T(n, d-1) + O(n) &"  se" n > 1 and d > 1
)
$

Dove $T(n, 1) = O(n)$ è il costo di costruzione dell'array associato (merge di due array ordinati), $2T(n/2, d)$ è il costo per costruire i sottoalberi sinistro e destro, $T(n, d-1)$ è il costo per costruire l'albero associato, e $O(n)$ è il costo per trovare il mediano e partizionare i punti.

La complessità di costruzione del layered range tree rimane $O(n log^(d-1) n)$ e lo stesso vale per lo spazio occupato.

== Query in un Layered Range Tree

L'algoritmo per eseguire una query ortogonale in un layered range tree è simile a quello del range tree fino alla penultima dimensione. Quando la si raggiunge, viene trovato il nodo di split e si esegue una ricerca binaria nell'array associato per trovare il primo elemento che è maggiore o uguale a $x_d$, quindi proseguendo nella ricerca nell'albero si seguono contemporaneamente i puntatori nell'array associato ad ogni nodo. Nel momento in cui si visita un nodo il cui sottoalbero è completamente contenuto nella regione di query, si riportano tutti gli elementi nell'array associato a partire dalla posizione puntata dal puntatore fino al più grande elemento che è minore o uguale a $x'_d$.


#figure(
  kind: "algorithm",
  supplement: [Algoritmo],

  pseudocode-list(booktabs: true, stroke:none, numbered-title: [Query in un Layered Range Tree])[
  + *function* LRT-RANGE-SEARCH(v, range, results):
    +  *if* v _è una foglia_:
      + *if* v.point _è in_ range:
        + aggiungi v.point a results
    + *else*:
      + _trova il nodo di split_ v_split

      + \// Segui il percorso verso x_n
      + v = v_split.left
      + *if* v.dim == d - 2:
        + pos = _ricerca binaria in_ v_split.assoc _per_ range.min[d-1]
      + *while* v _non è una foglia_:
        + dim = v.dim
        + *if* range.min[dim] < v.split:
          + \// Il sottoalbero destro è completamente contenuto nella regione di query
          + *if* dim == d - 2:
            + *while* pos < |v.right.assoc| *and* v.right.assoc[pos].point.coord[d-1] <= range.max[d-1]:
              + aggiungi v.right.assoc[pos].point a results
              + pos += 1
          + *else*:
            + LRT-RANGE-SEARCH(v.right.assoc, range, results)
          + v = v.left
          + *if* v.dim == d - 2:
            + pos = v.assoc[pos].left
        + *else*:
          + v = v.right
          + *if* v.dim == d - 2:
            + pos = v.assoc[pos].right
        + \// v è una foglia
        + *if* v.point _è in_ range:
            + aggiungi v.point a results
        + \
        + \// Segui il percorso verso x'\_n
        + \// ...
    ]
)

La ricerca binaria nell'array associato richiede $O(log n)$ tempo e viene eseguita solamente una volta, questa quindi non influisce nel calcolo della complessità in quanto la visita dei nodi occupa comunque O(log n) tempo. La ricerca nella struttura associata (della penultima dimensione) invece richiede solo il tempo necessario a riportare i punti trovati, eliminando così la ricerca all'interno dell'ultimo albero come avveniva nel range tree. Tenendo conto di questo la complessità della query diventa quindi:

$O(log^(d-1) n + k)$.

= Riassunto dei costi delle strutture dati

#figure(
  table(
    columns: 4,
    [*Struttura dati*], [*Costruzione*], [*Spazio*], [*Query*],
    [KD-Tree], [$O(n log n)$], [$O(n)$], [$O(n^(1-1/d) + k)$], 
    [Range Tree], [$O(n log^(d-1) n)$], [$O(n log^(d-1) n)$], [$O(log^d n + k)$], 
    [Layered Range Tree], [$O(n log^(d-1) n)$], [$O(n log^(d-1) n)$], [$O(log^(d-1) n + k)$], 
  )
)

= Implementazione (algoritmi.h)
Nella libreria fornita, le strutture dati sono implementate con classi C++ generiche per dimensione $D$. Per le primitive geometriche si utilizza CGAL (Computational Geometry Algorithms Library). Le classi che implementano le strutture dati sono:

`KDTree`: definisce tipi `InternalNode` e `LeafNode`. In un nodo interno si memorizzano `split_dimension` e `split_coordinate`, e `left`/`right` come `unique_ptr<Node>`. Viene usato `nth_element` su un vettore di punti per trovare il mediano sulla dimensione data. Le foglie (`LeafNode`) contengono direttamente un oggetto `Point<D>`. La funzione `range_search` controlla se un nodo è foglia (segna il punto se è nel range) oppure, in un nodo interno, decide se visitare i figli in base ai valori min/max della query lungo `split_dimension`.

`RangeTree`: definisce anch'esso `InternalNode` e `LeafNode`, dove ogni `InternalNode` ha `split_value`, `left` e `right`, più un puntatore `assoc` che fa riferimento ad un altro `RangeTree<D>` (per la dimensione successiva).

`LayeredRangeTree`: simile al `RangeTree`, ma l'associato (`assoc`) è di tipo generico `AssociatedStructure<D>`. Quando `dim+2 < D`, `assoc` è un'altro `LayeredRangeTree`; quando `dim+2 == D`, `assoc` è un `CascadeArray<D>`. Un cascade array è un vettore di elementi che memorizzano un punto e due indici opzionali `left` e `right` agli elementi corrispondenti nei figli.

Inoltre è stata implementata una visualizzazione interattiva dei kd-tree in due dimensioni usando `Qt5`.

#figure(
  image("imgs/visualization.png", width: 90%),
  caption: [Visualizzazione interattiva creata con Qt5 e CGAL di un kd-tree in 2D.]
)

== Benchmark

Per valutare le prestazioni delle strutture dati implementate, sono stati condotti benchmark su insiemi di punti bidimensionali uniformemente distribuiti di varie dimensioni (da 10.000 a 1.000.000 di punti). Per ogni dimensione del dataset, sono stati misurati:
- Tempo di costruzione (preprocessing)
- Utilizzo di memoria
- Tempo di query, come media su più query casuali con un numero costante di punti attesi nell'output

E' stato scelto di eseguire i benchmark con query che riportano un numero costante di punti per valutare meglio le differenze tra gli algoritmi, in quanto il tempo di query dipende linearmente dal numero di punti riportati.

Il kd-tree ha il tempo di costruzione minore e utilizza significativamente meno memoria, ma il tempo di query è notevolmente più alto degli altri algoritmi. Questo non è detto che sia il risultato della complessità peggiore dell'algoritmo, infatti data la piccola ampiezza del range delle query il termine $O(sqrt(n))$ occupa minor tempo rispetto al termine $O(log n)$ per la ricerca all'interno dell'albero.

Il Layered range tree, rispetto al range tree:
- utilizza a meno memoria, in quanto la struttura associata nell'ultima dimensione non memorizza i valori di split
- Ha un tempo di costruzione significativamente più basso, probabilmente dovuto al fatto che la costruzione dell'array associato è più veloce della costruzione di un albero di ricerca bilanciato
- Ha un tempo di query più basso, in quanto non deve eseguire ricerche binarie nell'ultima dimensione ed il impiegato nel riportare tutti i punti di un sottoalbero è maggiore rispetto al tempo impiegato nel riportare i punti in un array disposti in modo contiguo in memoria.

#figure(
  image("imgs/construction_time.png", width: 90%),
  caption: [Tempo di costruzione delle strutture dati al variare del numero di punti in scala log-log.])

#figure(
  image("imgs/memory_usage.png", width: 90%),
  caption: [Utilizzo di memoria delle strutture dati al variare del numero di punti in scala log-log.])

#figure(
  image("imgs/query_time.png", width: 90%),
  caption: [Tempo di query delle strutture dati al variare del numero di punti.])



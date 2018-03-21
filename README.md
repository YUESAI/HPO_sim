# HPO_sim

## Algorithm summary

  My algorithm is divided into two phases. The first phase is the division of the functional modules, as shown in the following figure:


![image](https://github.com/YUESAI/HPO_sim/raw/master/HPO_sim/Screenshots/1.jpg)

First, the module dictionary consisting of 349 function modules was extracted from the 13460 genes and 141,296 associated PPI networks obtained through functional modules. Among them, the function module extracts the PPI network-based Multi-tag Label Propagation Algorithm (Speaker-Listener Label Propagation Algorithm, SLPA). SLPA is an improved label propagation algorithm based on LPA (Label Propagation Algorithm, LPA).

LPA is a graph-based unsupervised learning method that uses topology information to assign labels to nodes. After the label is propagated, the nodes identified by the same label are divided into the same cluster.

But in our problem, a gene may belong to multiple functional modules, so I used the SLPA algorithm. Compared with LPA, the biggest feature of SLPA is that it records the historical tag sequence of each node during the refresh iteration, and when the iteration stops, the frequency of each (different) tag in each node's historical tag sequence is Statistics, according to a given threshold, filter out those tags that have a low frequency of occurrence, and the rest is the tag of the node (usually there are multiple).

***
The second stage of the algorithm is similarity calculation, as shown in the following figure:


![image](https://github.com/YUESAI/HPO_sim/raw/master/HPO_sim/Screenshots/2.jpg)

For example, I want to calculate the similarity between phenotype **A** and phenotype **B**. First, map their corresponding gene sets **Ga** and **Gb** through the module dictionary obtained in stage 1 to their respective module sets **Ma** and **Mb**, and then calculate the functional module set. The similarity between **Ma** and **Mb** gives similarity between phenotypes **A** and **B**.

When calculating the similarity of functional module collections, I used three indicators to calculate the similarity:


![image](https://github.com/YUESAI/HPO_sim/raw/master/HPO_sim/Screenshots/3.jpg)

## Key code summary

```
def get_mods(T, gene_set):

    mods = set()

    cnt = 0

    for i in gene_set:

        try:

            mods.update(T[i])

        except:

            mods.add(-1)

            cnt = cnt + 1

    # 对超过50%的基因不在字典中的情况拒判

    if cnt >= len(gene_set) * 0.5:

        raise '未能识别集合'

    return mods

 

def jac_sim(mods_1, mods_2):

    # Jaccard Index

    return len(mods_1 & mods_2) / len(mods_1 | mods_2)

 

def BP_sim(mods_1, mods_2):

    # 交集的大小

    l_and = len(mods_1 & mods_2)

    # Bin Prediction Accuracy

    return (l_and / len(mods_1) +  l_and / len(mods_2)) / 2

 

def aff_sim(mods_1, mods_2):

    # Affinity Score

return np.square(len(mods_1 & mods_2)) / len(mods_1) / len(mods_2)
```

## Analysis of results
The accuracy of each indicator on each test set is shown in the following table:

![image](https://github.com/YUESAI/HPO_sim/raw/master/HPO_sim/Screenshots/4.jpg)
## Advantages and disadvantages of the methods used
### Advantages
* We measure the phenotypic similarity using functional module set similarity rather than gene set similarity.
  * This is because in biological systems, genetically encoded proteins perform their functions in a modular manner through interactions while several papers (Bader and Hogue, 2003; Bu et al., 2003; Spirin and Mirny, 2003) also show that in networks and themselves Modules that are closely linked and less relevant to other parts generally correspond to meaningful biological units, such as protein complexes and functional modules.
  * Avoiding multiple genes with similar functions has too much influence on similarity, and the more specific genes are less affected by the lower number.
  * The collection of functional modules greatly reduces the size of the gene set and the efficiency of the algorithm is high.
* For the division of functional modules, the classes of the traditional community discovery algorithm LPA cannot be overlapped. The problem solved by SLPA is more in line with the project.


### Disadvantages
* The algorithm extracted by the function module still has a lot of room for improvement.
* Did not consider the difference in the influence of each module on the phenotype.

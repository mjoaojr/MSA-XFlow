from Bio import Phylo

tree = Phylo.read("BOX017.dnd", "newick")
print(tree)
Phylo.draw_ascii(tree)

print(f"Existem {tree.count_terminals()} sequências")
#print(f"Altura: {tree.depths(unit_branch_lengths=True).values()}")

def pre_ordem (tree, h):
    if (tree.is_preterminal()):
#        print (tree)
#        for i in tree.clades:
#            print (i)
        return h+1

#    print (tree)
    hmax = h
    for i in tree.clades:
        x = pre_ordem (i, h+1)
        hmax = max (hmax, x)

    return hmax

print (f"Altura: {pre_ordem (tree.clade, 1)}")

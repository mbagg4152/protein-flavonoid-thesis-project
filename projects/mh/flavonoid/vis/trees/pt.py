from Bio import Phylo as phylo
from Bio.Phylo import PhyloXMLIO as pxio
import pylab
fname = 'tree2.html'


# tree = phylo.read(fname, 'phyloxml')
# phx = pxio.read(fname)


def main():
    proc()


def proc():
    tree = phylo.read('agi-tree.xml', 'phyloxml')
    phylo.draw_graphviz(tree)
    pylab.show()

if __name__ == '__main__':
    main()

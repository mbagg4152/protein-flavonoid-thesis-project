from Bio import Phylo as phylo
from Bio.Phylo import PhyloXMLIO as pxio

fname = 'tree.html'
tree = phylo.read(fname, 'phyloxml')
phx = pxio.read(fname)


def main():
    proc()


def proc():
    phylo.draw_ascii(tree)


if __name__ == '__main__':
    main()

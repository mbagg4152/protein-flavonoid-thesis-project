from bioservices.kegg import KEGG
from bioservices import PDB
import os.path
import urllib.request, urllib.error, urllib.parse
import re

RE_SEQ = r'(&gt;.*</pre></div>)'
sep = os.sep
cwd = os.getcwd().replace('lib', '').replace('//', '/').replace('\\\\', '\\')
tst_dir = cwd + 'misc_output' + sep
kdb = KEGG()
pdata = PDB()

try: os.mkdir(tst_dir)
except FileExistsError: pass


def main():
    ans = True
    print('This is for testing kegg output format. output files wil be in misc_output' + sep)
    while ans:
        print(
            '\n1. Organism IDs' +
            '\n2. Output for cam00943' +
            '\n3. Get pathways for lja ' +
            '\n4. Look at dosa:Os08g0277200' +
            '\n5. Get 6DHL in PDB' +
            '\n6. Get 6DHL in XML' +
            '\n7. Get 6DHL in CIF' +
            '\n8. Get current PDB IDs' +
            '\n9. Get seq from https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+crb:17877590'
            '\n0. Exit'
            )
        ans = int(input('\nSelection: '))
        print('')
        if ans == 1:
            kdb.list('organisms')
            write_out('all-org-ids.txt', str(kdb.organismIds))
        elif ans == 2:
            raw = kdb.get('cam00941')
            raw_dict = kdb.parse(raw)
            gene = raw_dict.get('GENE', 'No data').values()
            compound = raw_dict.get('COMPOUND', 'No data').values()
            print('type of gene is ' + str(type(gene)) + ' type of compound is ' + str(type(compound)))
            write_out('pathway.txt', str(raw_dict))
            write_out('gene.txt', str(gene))
            write_out('compound.txt', str(compound))
        elif ans == 3:
            kdb.organism = 'lja'
            write_out('single-org.txt', str(kdb.pathwayIds))
        elif ans == 4:
            raw = kdb.get('dosa:Os08g0277200')
            raw_dict = kdb.parse(raw)
            gene = raw_dict.get('NTSEQ', 'No data')
            print(str(gene))
        elif ans == 5:
            res = pdata.get_file('6DHL', 'pdb')
            write_out('6dhl.pdb', res)
        elif ans == 6:
            url = "https://files.rcsb.org/view/6DHL.xml"
            urllib.request.urlretrieve(url, tst_dir + '6DHL.xml')
        elif ans == 7:
            res = pdata.get_file('6DHL', 'cif')
            write_out('6dhl.cif', res)
        elif ans == 8:
            res = pdata.get_current_ids()
            write_out('current-ids.txt', str(res))
        elif ans == 9:
            try: urllib.request.urlretrieve('https://www.kegg.jp/dbget-bin/www_bget?-f+-n+n+crb:17877590',
                                            tst_dir + 'seq.txt')
            except urllib.error.HTTPError or urllib.error.URLError as e:
                print('err getting seq')
                return

            try:
                file = open(tst_dir + 'seq.txt', 'r')
                info = file.read().replace('\\n', '').replace('\n', ' ')
                seq = re.findall(RE_SEQ, info)

                try:
                    out = seq[0].replace('&gt;', '>').replace('</pre></div>', '')
                    dna = re.findall(r'([atcg]{8,})', out)
                    dnaf = ''.join(dna)
                    print(out)
                    print(dnaf)
                except IndexError: print('mustve not found what was needed')


            except FileNotFoundError: print('couldnt find file')


def write_out(name, contents):
    try:
        wf = open(tst_dir + name, 'x')
        wf.close()
    except FileExistsError: print('')
    wf = open(tst_dir + name, 'w')
    wf.write(contents)
    wf.close()
    print('Wrote to: ' + name)


if __name__ == '__main__':
    main()

from bioservices.kegg import KEGG
import os.path

sep = os.sep
cwd = os.getcwd()
tst_dir = cwd + sep + 'misc_output' + sep
kdb = KEGG()
fn_org_ids = 'allOrganismIDs.txt'
fn_org_pathway = 'orgPathGetOutput.txt'
fn_single_org_paths = 'pathsForSingleOrg.txt'
fn_gene = 'geneDataOnly.txt'
fn_compound = 'compoundDataOnly.txt'


def main():
    ans = True
    print('This is for testing kegg output format. output files wil be in misc_output' + sep)
    while ans:
        print(
            '\n1. Organism IDs' +
            '\n2. Output for cam00943' +
            '\n3. Get pathways for lja ' +
            '\n4. Look at ats:LOC109771219' +
            '\n0. Exit'
        )
        ans = int(input('\nSelection: '))
        print('')
        if ans == 1:
            kdb.list('organisms')
            write_out(fn_org_ids, str(kdb.organismIds))
        elif ans == 2:
            raw = kdb.get('cam00941')
            raw_dict = kdb.parse(raw)
            gene = raw_dict.get('GENE', 'No data').values()
            compound = raw_dict.get('COMPOUND', 'No data').values()
            print('type of gene is ' + str(type(gene)) + ' type of compound is ' + str(type(compound)))
            write_out(fn_org_pathway, str(raw_dict))
            write_out(fn_gene, str(gene))
            write_out(fn_compound, str(compound))
        elif ans == 3:
            kdb.organism = 'lja'
            write_out(fn_single_org_paths, str(kdb.pathwayIds))
        elif ans == 4:
            raw = kdb.get('ats:LOC109771219')
            raw_dict = kdb.parse(raw)
            gene = raw_dict.get('NTSEQ', 'No data')
            print(str(type(gene)))

def write_out(name, contents):
    wf = open(tst_dir + name, 'w')
    wf.write(contents)
    wf.close()
    print('Wrote to: ' + name)


if __name__ == '__main__':
    main()

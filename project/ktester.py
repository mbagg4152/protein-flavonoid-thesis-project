from bioservices.kegg import KEGG
import os.path

sep = os.sep
cwd = os.getcwd()
tst_dir = cwd + sep + 'test-output' + sep
kdb = KEGG()
fn_org_ids = 'allOrganismIDs.txt'
fn_org_pathway = 'orgPathGetOutput.txt'
fn_single_org_paths = 'pathsForSingleOrg.txt'


def main():
    ans = True
    print('This is for testing kegg output format. output files wil be in test-output' + sep)
    while ans:
        print(
            '\n1. Organism IDs' +
            '\n2. Output for cam00940' +
            '\n3. Get pathways for lja ' +
            '\n4.' +
            '\n0. Exit'
        )
        ans = int(input('\nSelection: '))
        print('')
        if ans == 1:
            kdb.list('organisms')
            write_out(fn_org_ids, str(kdb.organismIds))
        elif ans == 2:
            write_out(fn_org_pathway, kdb.get('cam00940'))
        elif ans == 3:
            kdb.organism = 'lja'
            write_out(fn_single_org_paths, str(kdb.pathwayIds))


def write_out(fname, contents):
    # full_name = os.path.join(tst_dir, fname)
    wf = open(tst_dir + fname, 'w')
    wf.write(contents)
    wf.close()
    print('Wrote to: ' + fname)


if __name__ == '__main__':
    main()

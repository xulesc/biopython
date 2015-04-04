# Copyright (C) 2015, Anuj Sharma (anuj.sharma80@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
This module provides an implementation for calculating similarity between two
protein structures using TM-align (tmscore). Description and theoretical basis
for its applicability to protein structure comparison can be found in

Y. Zhang, J. Skolnick, "TM-align: A protein structure alignment algorithm based
on TM-score", Nucleic Acids Research, 33: 2302-2309 (2005)
"""

import os.path as opath
from subprocess import check_output as chkout
from Bio.PDB.PDBList import PDBList


class TMalign:

    """
    TMalign is a method for pairwise comparison of protein structure. The
        method identifies alignments between protein pairs combining the TM-score
        rotation matrix and Dynamic Programming (DP). Three types of initial
        alignments are generated from the secondary structures of the protein. The
        initial alignments are submitted to a heuristic iterative algorithm for
        refining the alignment.

    Reference:
        Y. Zhang, J. Skolnick, "TM-align: A protein structure alignment algorithm
        based on TM-score", Nucleic Acids Research, 33: 2302-2309 (2005)
        """

        def __init__(self, prog='tmalign'):
        """
        Initialize the TMalign object. Default values for the TMalign program
		is 'tmalign' in the current working directory. The location of the
		program can be overridden to point to other location where the program
		may have been downloaded.
        """
            self._prog = './' + prog
            if not opath.isfile(self._prog):
                raise IOError('cannot open', prog)
            self._pdbl = PDBList()
            self._alternator = self._alternate()
            self._read_aln = False
            self._sequences = ['', '']
            # setup parsers for the text level interface with TMalign
            glen = lambda x: x.split()[-2]
            gname = lambda x: '\'' + x.split()[-1] + '\''
            gtm = lambda x: x.split()[1]
            gdet = lambda x: str(x.replace('Aligned length',
                                           'Alignedlength').replace(',', '').split()[1:6:2])
            self._txt2dat = {'Name of Chain_1': (gname, 'self._prot1_name'),
                             'Name of Chain_2': (gname, 'self._prot2_name'),
                             'Length of Chain_1': (glen, 'self._prot1_len'),
                             'Length of Chain_2': (glen, 'self._prot2_len'),
                             'Chain_1)': (gtm, 'self._tm1'), 'Chain_2)': (gtm, 'self._tm2'),
                             ' length=': (gdet, '(self._allen, self._rmsd, self._seqid)'),
                             '(":"': (lambda x: 'True', 'self._read_aln')}

        def _alternate(self):
            """
            Generator that alternates return values of 0 and 1
            """
            while True:
                yield 0
                yield 1

        def _stack_alignments(self, line):
            """
            Incrementally builds the alignment string for the two structures. Odd
            iterations result in increamenting the sequence for the first protein.
            Even iterations result in incrementing the sequnce for the second
            protein.
            """
            if len(line) == 0 or line[0] in [' ', '.', ':']:
                return
            self._sequences[self._alternator.next()] += line

        def _parse_tmalign_output(self):
            """
            This is the text level interface with TMalign. The output of the program
            is parsed and relevant values of the object are populated.
            """
            for oline in self._output.split('\n'):
                for k, v in self._txt2dat.iteritems():
                    if self._read_aln:
                        self._stack_alignments(oline)
                        break
                    elif k in oline:
                        exec(v[1] + ' = ' + v[0](oline))
                        break
            self._read_aln = False

        def _run_tmalign(self, p1file, p2file):
            """
            Tests for availability of the two pdb files and calls the external
            TMalign program along with the two pdb files as parameters.
            """
            if not opath.isfile(p1file):
                raise IOError('cannot open', p1file)
            if not opath.isfile(p2file):
                raise IOError('cannot open', p2file)
            self._sequences = ['', '']
            self._output = chkout([self._prog, p1file, p2file])
            self._parse_tmalign_output()

        def frun(self, p1file, p2file):
            """
            Calls the external TMalign program along with the two pdb files as
            parameters.
            """
            self._run_tmalign(p1file, p2file)

        def run(self, pdb1='1AA9', pdb2='1ASH', odir='/tmp'):
            """
            Calls the external TMalign program along with the two pdb files as
            parameters. The pdb files are created by downloading the domains from
            pdb. The files are downloaded to the directory specified by the odir
            argument.
            """
            p1file = self._pdbl.retrieve_pdb_file(pdb1, pdir=odir)
            p2file = self._pdbl.retrieve_pdb_file(pdb2, pdir=odir)
            self.frun(p1file, p2file)

        def getTMScore1(self):
            return self._tm1

        def getTMScore2(self):
            return self._tm2

        def getRMSD(self):
            return self._rmsd

        def getAlignLength(self):
            return self._allen

        def getSequenceID(self):
            return self._seqid

        def getProt1Len(self):
            return self._prot1_len

        def getProt2Len(self):
            return self._prot2_len

        def getAlignment(self):
            return self._sequences

        def usageTMalign(self):
            return chkout([self._prog])

if __name__ == '__main__':
    tmalign = TMalign()
    tmalign.run()
    print tmalign.getTMScore1()
    print tmalign.getTMScore2()
    print tmalign.getRMSD()
    print tmalign.getAlignment()

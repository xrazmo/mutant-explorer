import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class gemucator(object):

    '''The gemucator class (short for GenBank Mutation Locator) is designed to tell you the location of a mutation
    in the given GenBank reference genome.

    Notes:
        * Mutations are specified in the form gene_mutation e.g. rpoB_S450L. Three forms are accepted: SNP, PROMOTER or INDEL.
        * Amino acids are specified throughout in UPPERCASE and nucleotides in lower case.
        * A SNP mutation specifies the reference and target amino acids and the residue position in the coding sequence e.g. katG_S315T.
        * A PROMOTER mutation requires the reference and target nucleotides to be specified, as well as the location relative to the start codon e.g. inhA_c-15t
        * An INDEL mutation is less specific and requires the location, whether an insertion (ins) or deletion (del) is occuring, and the number of bases (or * as a wildcard for 'any number') e.g. ahpC_66_del_1
        * The code is extremely defensive and if the mutation nomenclature includes the reference amino acid or base, this will be checked against the provided genbank file.
    '''

    def __init__(self, genbank_file):
        '''
        Instantiate an instance of the class by loading the specified GenBank file

        Args:
            genbank_file (str) path to the required GenBank file. By default version 3 of the H37rV genome for M. tuberculosis is loaded.
        '''

        # use BioPython to load the genbank file
        self.genome=SeqIO.read(genbank_file,'genbank')

        # remember the name of the genbank file for use in assert statements later
        self.genbank_file=genbank_file

        self.version=self.genome.annotations['accessions'][0]+"."+str(self.genome.annotations['sequence_version'])

    def locate_mutation(self,mutation,nucleotide_mutation=False):

        '''
        Given a mutation, return the numeric location in the GenBank file.

        Args:
            mutation (str) e.g. katG_S315T, katG_c-15t, katG_200_ins_3
            nucleotide_mutation (bool). Set to True if this is an RNA encoding gene (rather than encoding amino acids)

        Returns:
            A tuple.

            If the mutation does not validate e.g. katG_F315T, then this has the form (False,'S'). The second element is the base or amino acid at this position in the Genbank file,
            If it does validate then for a nucleotide or amino acid mutation the tuple is (True,('c',1673424)) or (True,('agc', [2155169, 2155168, 2155167]))
        '''

        result=self._analyse_mutation(mutation,nucleotide_mutation=nucleotide_mutation)

        return(result)
        # if result[0]:
        #     return(result[1])
        # else:
        #     return(False)

    def valid_mutation(self,mutation):

        '''
        Simply checks to see if the specified mutation validates against the supplied reference genbank file.

        Args:
            mutation (str) e.g. katG_S315T, katG_c-15t, katG_200_ins_3
            nucleotide_mutation (bool). Set to True if this is an RNA encoding gene (rather than encoding amino acids)

        Returns:
            True/False
        '''

        # find out what the type of gene it is
        gene_name=mutation.split("_")[0]

        # find out if it is a GENE, LOCUS or RNA
        gene_type=self.gene_type(gene_name)

        # if it is a gene encoding rRNA, then tell the method!
        if gene_type=="RNA":
            result=self._analyse_mutation(mutation,nucleotide_mutation=True)
        else:
            result=self._analyse_mutation(mutation,nucleotide_mutation=False)

        if result[0]:
            return(True)
        else:
            return(False)

    def _parse_mutation(self,mutation,nucleotide_mutation):
        '''
        Parse the mutation and return a collection of variables and Booleans.

        Args:
            mutation (str) e.g. katG_S315T, katG_c-15t, katG_200_ins_3

        Returns:
            gene_name (str): the name of the gene
            before (str): the amino acid or base in the reference genomes
            position (int or str): the numerical position of either the amino acid or base
            after (str): the mutated amino acid or base
            wildcard (bool): whether the mutation applies to all positions or not
            promoter (bool): whether the mutation applies to the promoter of a gene
            mutation_type (str): either PROMOTER, SNP or INDEL
        '''

        # first, parse the mutation
        cols=mutation.split("_")

        before=None
        after=None
        wildcard=False
        promoter=False

        # check the mutation at least comprises the expected number of sections
        assert len(cols) in [2,3,4], "mutation "+mutation+" not in correct format!"

        # the gene/locus name should always be the first component
        gene_name=cols[0]

        assert self.valid_gene(gene_name), "gene not found in GenBank file! "+gene_name

        # determine if this is a CDS or PROMOTER SNP mutation
        if len(cols)==2:

            if '*' in cols[1]:

                # there can be no 'before' amino acid if there is a wildcard at position
                wildcard=True

                if cols[1][0]=='-':
                    promoter=True
                    position=str(cols[1][1:-1])
                else:
                    promoter=False
                    position=str(cols[1][0:-1])

                # all the positions should be a wildcard otherwise something is wrong..
                assert position=='*', mutation+' has a * but not formatted like a wildcard'
            else:

                wildcard=False

                before=cols[1][0]
                position=int(cols[1][1:-1])

                if position<0:
                    promoter=True
                else:
                    promoter=False

            # they all have the after amino acid in the same place
            after=cols[1][-1]

            # if it is a promoter mutation
            if promoter:

                if not wildcard:
                    assert before in ['c','t','g','a'], before+" is not a nucleotide!"

                assert after in ['c','t','g','a','?','z'], after+" is not a nucleotide!"

                mutation_type="PROMOTER"

            # ..otherwise it is an amino acid SNP
            else:

                if nucleotide_mutation:

                    if not wildcard:
                        assert before in ['c','t','g','a'], before+" is not a nucleotide!"

                    assert after in ['c','t','g','a','?','z'], after+" is not a nucleotide!"

                else:

                    if not wildcard:
                        assert before in ["!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], before+" is not an amino acid!"

                    assert after in ['=','?',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"

                mutation_type="SNP"

        # otherwise it must be an INDEL, which is always nucleotide based
        else:

            # deal with wildcards in the position
            if cols[1] in ["*","-*"]:

                wildcard=True

            else:

                wildcard=False

                position=int(cols[1])

                # be defensive here also!
                assert cols[2] in ["ins","del","indel","fs"], "INDEL must be on the format rpoB_1300_ins_1 i.e. the third element must be ins or del, not "+cols[2]

                if len(cols)==4:
                    if cols[3].isnumeric():
                        assert int(cols[3])>0, "number of nucleotides inserted or deleted must be >0"
                    else:
                        assert bool(re.match('^[catg]+$', cols[3])), "INDEL contains bases other than a,t,c,g"

                # only then allow this to be an INDEL!
                mutation_type="INDEL"

        return(gene_name,before,position,after,wildcard,promoter,mutation_type)


    def _analyse_mutation(self,mutation,nucleotide_mutation=False):

        '''Given a specified mutation, return the location(s) and nucleotide(s) in the reference genome.

        Args:
            mutation (str): a SNP, PROMOTER or INDEL mutation in the format described in the class docstring.
            nucleotide_mutation (bool): if True, forces the code to consider this a simple nucleotide mutation (i.e. not an amino acid). Useful for RNA genes like rrs.

        Returns:
            result (boolean): whether successful or not
            (bases, base_positions): a tuple of the bases and their locations. Hence an amino acid mutation would have 3 entries.
        '''

        # first parse the supplied mutation
        (gene_name,before,position,after,wildcard,promoter,mutation_type)=self._parse_mutation(mutation,nucleotide_mutation)

        # now look up the gene in the genbank reference
        (result,record)=self._fetch_genbank_record(gene_name)


        base_positions=[]

        # if
        if not result:
            return(False,(None,None))

        else:

            # if we are here the gene exists, but there is a wildcard so we can't check any further
            if wildcard:
                return(True,(None,None))

            else:

                # the start and end positions in the reference genome
                start=record.location.start.position
                end=record.location.end.position

                # and finally whether the strand (which if -1 means reverse complement)
                strand=record.location.strand.real

                # retrieve the coding nucleotides for this gene
                coding_nucleotides=self.genome[start:end]

                if mutation_type=="SNP" and not nucleotide_mutation:

                    if strand==1:

                        base_position=start-1+((3*position)-2)

                        bases=self.genome[base_position:base_position+3].seq

                        if before=="!":
                            if "*"!=bases.translate():
                                return(False,str(bases.translate()))
                        else:
                            if before!=bases.translate():
                                return(False,str(bases.translate()))

                        bases=str(bases).lower()

                        for i in bases:
                            base_positions.append(base_position)
                            base_position+=1

                    else:

                        base_position=end+1-((3*position)-2)

                        bases=self.genome[base_position-3:base_position].reverse_complement().seq

                        if before=="!":
                            if "*"!=bases.translate():
                                return(False,str(bases.translate()))
                        else:
                            if before!=bases.translate():
                                return(False,str(bases.translate()))

                        bases=str(bases).lower()

                        for i in bases:
                            base_positions.append(base_position)
                            base_position+=-1

                elif (mutation_type in ["PROMOTER","INDEL"]) or nucleotide_mutation:

                    if strand==1:

                        if mutation_type=="PROMOTER":
                            base_position=start+position
                        else:
                            base_position=start+position-1

                        bases=self.genome[base_position:base_position+1].seq.lower()

                    else:

                        if mutation_type=="PROMOTER":
                            base_position=end-position-1
                        else:
                            base_position=end-1-position

                        bases=self.genome[base_position:base_position+1].reverse_complement().seq.lower()

                    bases=str(bases).lower()

                    # as the reference base is specified for a promoter mutation, be defensive and check it agrees with the genbank file
                    if mutation_type=="PROMOTER":
                        if bases!=before:
                            return(False, str(bases))

                    base_positions.append(base_position)

                else:

                    raise (False,"are you sure mutation "+mutation+" is in this reference genome?")

        if not base_positions:
            raise Exception("are you sure mutation "+mutation+" is in this reference genome?")

        return(True,(bases,base_positions))

    def valid_gene(self,gene_name):

        '''
        Simply checks to see if the specified gene exists in the supplied reference genbank file by searching all the gene and locus tags.

        Args:
            gene_name (str) e.g. katG

        Returns:
            True/False
        '''

        result=self._fetch_genbank_record(gene_name)

        if result[0] is False:
            return(False)
        else:
            return(True)

    def gene_type(self,gene_name):

        '''
        Returns GENE or LOCUS as the gene type. Cannot distinguish RNA genes at present.

        Args:
            gene_name (str) e.g. "katG"

        Returns:
            GENE, LOCUS, RNA or None
        '''

        result=self._fetch_genbank_record(gene_name)

        if result[0] is False:
            return(None)
        else:
            if result[1].type=="rRNA":
                return("RNA")
            else:
                return(result[0])

    def _fetch_genbank_record(self,gene_name):

        '''
        Internal method that feeds both gene_type() and valid_gene()

        Args:
            gene_name (str) e.g. katG

        Returns:
            Tuple with first element being GENE,LOCUS or None and the second being the record in the GenBank file
        '''

        assert gene_name is not None, "you have to specify a gene_name!"

        for record in self.genome.features:

            # check that the record is a Coding Sequence
            if record.type in ['CDS','rRNA']:

                found_record=False

                # if it is a gene
                if 'gene' in record.qualifiers.keys():
                    if gene_name in record.qualifiers['gene']:
                        return('GENE',record)

                elif 'locus_tag' in record.qualifiers.keys():
                    if gene_name in record.qualifiers['locus_tag']:
                        return('LOCUS',record)

        return(False,None)

    def identify_gene(self,location,alt_base=None,promoter_length=100):

        '''
        What is encoded at the supplied numeric position in the genome defined by the genbank file?

        Args:
            location (int)  numeric location in the genome
            promoter_length (int) the number of nucleotides upstream of the coding sequence to consider as the promoter. Default is 100.

        Returns:
            gene_name (str) the name of the name, if relevant, e.g. "katG"
            residue or base (str) e.g. "A" or "a"
            position (int)
        '''
       
        assert int(location)>0, "genomic position has to be a positive integer"

        assert int(location)<len(self.genome), "genome position has to be less than the size of the genome...!"

        residue,gene_name,ref_base,alt_residue,product,prot_id = None,None,None,None,None,None
        # now that we know what we are looking for, iterate through all the features in the genomes
        for record in self.genome.features:

            found_gene=False

            # check that the record is a Coding Sequence
            if record.type in ['CDS','rRNA']:

                # the start and end positions in the reference genome
                start=record.location.start.position
                end=record.location.end.position
                strand=record.location.strand.real

                # regular strand
                if strand==1:

                    # are we within the expected region?
                    if (start-promoter_length) < location < end:

                        found_gene=True
                        ref_base = self.genome[location].lower()
                        # retrieve the coding nucleotides for this gene
                        coding_nucleotides=self.genome[start:end]

                        # find out the gene_name
                        if 'gene' in record.qualifiers.keys():
                            gene_name=record.qualifiers['gene'][0]
                        elif 'locus_tag' in record.qualifiers.keys():
                            gene_name=record.qualifiers['locus_tag'][0]
                        else:
                            gene_name=None
                        
                        if 'product' in record.qualifiers.keys():
                            product=record.qualifiers['product'][0]
                        
                        if 'protein_id' in record.qualifiers.keys():
                            prot_id=record.qualifiers['protein_id'][0]
                        # if we are in the CDS
                        if location>start:

                            if self.gene_type(gene_name)!="RNA":
                                position=(location-start+3)//3
                                residue=coding_nucleotides.seq.translate()[position-1]

                                if alt_base != None:
                                    #replacing the base with the point mutation
                                    tmp_genome = list(str(self.genome.seq))
                                    tmp_genome[location] = alt_base
                                    coding_nucleotides = SeqRecord(Seq(''.join(tmp_genome[start:end])),id=self.genome.id)
                                    alt_residue = coding_nucleotides.seq.translate()[position-1]
                            else:
                                position=location-start+1
                                residue=self.genome[location].lower()
                                
                                

                        # otherwise we are a promoter
                        else:

                            position=location-start
                            residue=self.genome[location].lower()
                            

                elif strand==-1:

                    # are we within the expected region?
                    if start < location < end+promoter_length:

                        found_gene=True
                        ref_base = self.genome[location].lower()
                        # retrieve the coding nucleotides for this gene
                        coding_nucleotides=self.genome[start:end]

                        # find out the gene_name
                        if 'gene' in record.qualifiers.keys():
                            gene_name=record.qualifiers['gene'][0]
                        elif 'locus_tag' in record.qualifiers.keys():
                            gene_name=record.qualifiers['locus_tag'][0]
                        else:
                            gene_name=None

                        if 'product' in record.qualifiers.keys():
                            product=record.qualifiers['product'][0]
                        
                        if 'protein_id' in record.qualifiers.keys():
                            prot_id=record.qualifiers['protein_id'][0]

                        if location<end:

                            if self.gene_type(gene_name)!="RNA":
                                position=(end-location+3)//3
                                residue=coding_nucleotides.reverse_complement().seq.translate()[position-1]
                                
                                if alt_base != None:
                                    #replacing the base with the point mutation
                                    tmp_genome = list(str(self.genome.seq))
                                    tmp_genome[location] = alt_base
                                    coding_nucleotides = SeqRecord(Seq(''.join(tmp_genome[start:end])),id=self.genome.id)
                                    alt_residue=coding_nucleotides.reverse_complement().seq.translate()[position-1]
                            else:
                                position=end-location+1
                                residue=self.genome[location:location+1].reverse_complement().seq.lower()

                        # otherwise we are a promoter
                        else:

                            position=end-location
                            residue=self.genome[location:location+1].reverse_complement().seq.lower()


                if found_gene:
                    return{"gene_name":gene_name,
                            "product":product,
                            "protein_id":prot_id,
                            "pos_gene":position,
                            "ref_aa":residue,
                            "alt_aa":alt_residue,
                            "ref_nt2":ref_base,
                            "strand":strand}

        return {}
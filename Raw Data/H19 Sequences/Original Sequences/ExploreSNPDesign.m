
H19_seq = fastaread('AK145379.exons.fasta');
H19_seq2 = readfasta('AK145379.exons.fasta');

B6_SNP0 = 'ccagcccatggtgttcaagaaggc';
C7_SNP0 = 'ccagtccatggtgttcaagaaggc';

B6_SNP1 = 'actcaaagctatctccgggactcc';
C7_SNP1 = 'actccaagctatctccgggactcc';

B6_SNP2 = 'gtttacacactcgctgtatacattcatac';
C7_SNP2 = 'gtttgcacactcgctgtatacattcatac';

B6_SNP4 = 'tggacgacaggtgggtactgggg';
C7_SNP4 = 'tggatgacaggtgggtactgggg';


%%

SNPs = cell(9,1);

SNPs{1} = B6_SNP0;
SNPs{2} = C7_SNP0;
SNPs{3} = B6_SNP1;
SNPs{4} = C7_SNP1;
SNPs{5} = B6_SNP2;
SNPs{6} = C7_SNP2;
SNPs{7} = B6_SNP4;
SNPs{8} = C7_SNP4;
SNPs{9} = H19_seq2;

SNPs = revcomp(SNPs)

%%
SNPs_B6 = SNPs([1,3,5,7]);
SNPs_C7 = SNPs([2,4,6,8]);

SNPs_B6 = cell2struct(SNPs_B6);
SNPs_C7 = cell2struct(SNPs_C7);
%%

test = localalign(H19_seq2, 'AGCCTTCTTGAACACCATGGGCTGGC')

test.Alignment{1}
 
 %%
test = multialign(SNPs);
seqalignviewer(test)

%%

test = nwalign(H19_seq2, SNPs{1})



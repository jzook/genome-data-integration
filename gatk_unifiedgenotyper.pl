#read in user parameters
my $in = shift @ARGV;  #sorted bam file (without ".bam")
my $res=0;

#$res=0;
if ($res == 0) {
    my $samCmd = "java -jar -Xmx20g  ~/bioinfo/GenomeAnalysisTK-1.6-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T UnifiedGenotyper -I ${in}.bam -o ${in}.gatk.conf2.mbq10.snpindel.raw.vcf -stand_call_conf 2.0 -stand_emit_conf 2.0 -glm BOTH -mbq 10 -A AlleleBalance --num_threads 8";
    print "$samCmd\n";
    $res = system($samCmd);
}   

if ($res == 0) {
    my $samCmd = "java -jar -Xmx2g  ~/bioinfo/GenomeAnalysisTK-1.4-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T SelectVariants -V ${in}.gatk.conf2.mbq10.snpindel.raw.vcf -o ${in}.gatk.conf2.mbq10.snp.raw.vcf -selectType SNP &";
    print "$samCmd\n";
    $res = system($samCmd);
}   

if ($res == 0) {
    my $samCmd = "java -jar -Xmx2g  ~/bioinfo/GenomeAnalysisTK-1.4-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T LeftAlignVariants -V ${in}.gatk.conf2.mbq10.snpindel.raw.vcf -o ${in}.gatk.conf2.mbq10.snpindel.LAlign.raw.vcf";
    print "$samCmd\n";
    $res = system($samCmd);
}   

if ($res == 0) {
    my $samCmd = "java -jar -Xmx2g  ~/bioinfo/GenomeAnalysisTK-1.4-5/dist/GenomeAnalysisTK.jar -R /data/results/justin/GRCh37/human_g1k_v37.fasta -T SelectVariants -V ${in}.gatk.conf2.mbq10.snpindel.LAlign.raw.vcf -o ${in}.gatk.conf2.mbq10.indel.raw.vcf -selectType INDEL";
    print "$samCmd\n";
    $res = system($samCmd);
}   


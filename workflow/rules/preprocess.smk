# Copyright 2024 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


rule extract_1KG_biallelic_snps:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
    output:
        vcf = "results/processed_data/1KG/snps/1KG.chr{i}.biallelic.snps.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -m 2 -M 2 -v snps | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_1KG_anc_alleles:
    input:
        vcf = rules.extract_1KG_biallelic_snps.output.vcf,
    output:
        info = "results/processed_data/1KG/anc_alleles/1KG.chr{i}.anc.alleles.bed",
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%POS\\t%INFO/AA\\n" {input.vcf} | awk '$4!~/\\.\\|\\|\\|/' | sed 's/|//g' |sed 's/[atcg]/\\U&/' | awk '$4~/[ATCG]/' | awk '{{print $1"\\t"$2"\\t"$3+1"\\t"$4}}' > {output.info}
        """


rule merge_1KG_src:
    input:
        vcf = rules.extract_1KG_biallelic_snps.output.vcf, 
        nea_vcf = rules.download_nea_genome.output.vcf,
        den_vcf = rules.download_den_genome.output.vcf,
    output:
        vcf = "results/processed_data/1KG/merged/merged_1KG_{src}.chr{i}.vcf.gz",
    resources:
        time = 4320,
    shell:
        """
        if [ {wildcards.src} == "nea" ]; then
            bcftools merge {input.vcf} {input.nea_vcf} | bcftools view -v snps -m 2 -M 2 -i "INFO/AN=5010" | bgzip -c > {output.vcf}
        else
            bcftools merge {input.vcf} {input.den_vcf} | bcftools view -v snps -m 2 -M 2 -i "INFO/AN=5010" | bgzip -c > {output.vcf}
        fi
        tabix -p vcf {output.vcf}
        """


rule merge_1KG_nea_den:
    input:
        vcf1 = "results/processed_data/1KG/merged/merged_1KG_nea.chr{i}.vcf.gz", 
        vcf2 = rules.download_den_genome.output.vcf,
    output:
        vcf = "results/processed_data/1KG/merged/merged_1KG_nea_den.chr{i}.vcf.gz",
    resources:
        time = 4320,
    shell:
        """
        bcftools merge {input.vcf1} {input.vcf2} | bcftools view -v snps -m 2 -M 2 -i "INFO/AN=5012" | bgzip -c > {output.vcf}
        """


rule create_ref_tgt_samples:
    input:
        samples = rules.download_1KG_info.output.samples,
    output:
        ref = "results/processed_data/1KG/samples/1KG.afr_eas_ref.samples.txt",
        tgt = "results/processed_data/1KG/samples/1KG.eur_tgt.samples.txt",
    shell:
        """
        grep -w 'AFR' {input.samples} | grep -v ACB | grep -v ASW | awk '{{print "AFR+EAS\\t"$1}}' > {output.ref}
        grep -w 'EAS' {input.samples} | awk '{{print "AFR+EAS\\t"$1}}' >> {output.ref}
        grep -w 'EUR' {input.samples} | awk '{{print $3"\\t"$1}}' > {output.tgt}
        """


rule annotate_1KG:
    input:
        vcf = rules.merge_1KG_nea_den.output.vcf,
        avsnp150 = rules.download_annovar_db.output.avsnp150,
        dbnsfp42c = rules.download_annovar_db.output.dbnsfp42c,
    output:
        avinput = "results/annotated_data/1KG/merged_1KG_nea_den.chr{i}.avinput",
        txt = "results/annotated_data/1KG/merged_1KG_nea_den.chr{i}.hg19_multianno.txt",
    resources:
        cpus=8, mem_gb=32,
    params:
        output_prefix = "results/annotated_data/1KG/merged_1KG_nea_den.chr{i}",
    shell:
        """
        bcftools query -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\n" {input.vcf} > {output.avinput}
        resources/tools/annovar/table_annovar.pl {output.avinput} resources/tools/annovar/humandb/ -buildver hg19 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . --thread {resources.cpus}
        """
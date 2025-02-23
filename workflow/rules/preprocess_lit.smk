# Copyright 2025 Xin Huang
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


ruleorder: beagle_imputation_with_ref > extract_chromosome


rule plink_to_vcf:
    input:
        ped = rules.download_lit_genomes.output.ped,
        map = rules.download_lit_genomes.output.map,
    output:
        vcf = "results/processed_data/lit/lit.vcf.gz",
    shell:
        """
        plink --ped {input.ped} --map {input.map} --recode vcf --out results/processed_data/lit/lit
        bgzip -c results/processed_data/lit/lit.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_chromosome:
    input:
        vcf = rules.plink_to_vcf.output.vcf,
    output:
        vcf = "results/processed_data/lit/chromosomes/lit.chr{i}.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -r {wildcards.i} | bgzip -c > {output.vcf}
        """


rule align_vcf_to_reference:
    input:
        vcf = rules.extract_chromosome.output.vcf,
        reference = rules.download_1KG_info.output.ref,
    output:
        vcf = "results/processed_data/lit/chromosomes/lit.chr{i}.aligned.vcf.gz",
    script:
        "../scripts/align_reference.py"


rule tabix_index:
    input:
        vcf = rules.align_vcf_to_reference.output.vcf,
    output:
        vcf = "results/processed_data/lit/chromosomes/lit.chr{i}.aligned.vcf.gz.tbi",
    shell:
        """
        tabix -p vcf {input.vcf}
        """


rule extract_yoruba_eur_samples:
    input:
        samples = rules.download_1KG_info.output.samples,
    output:
        yri_samples = "results/processed_data/lit/yri_samples.txt",
        eur_samples = "results/processed_data/lit/european_samples.txt",
    shell:
        """
        grep "YRI" {input.samples} | awk '{{print $1}}' > {output.yri_samples}
        grep -E "CEU|TSI|FIN|GBR|IBS" {input.samples} | awk '{{print $1}}' > {output.eur_samples}
        """


rule filter_for_eur:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
        samples = rules.extract_yoruba_eur_samples.output.eur_samples,
    output:
        vcf = "results/processed_data/lit/panels/chr{i}.EUR.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -S {input.samples} -v snps -m 2 -M 2 | bcftools norm -d none | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule filter_for_yri:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
        samples = rules.extract_yoruba_eur_samples.output.yri_samples,
    output:
        vcf = "results/processed_data/lit/panels/chr{i}.YRI.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -S {input.samples} -v snps -m 2 -M 2 | bcftools norm -d none | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule beagle_imputation_with_ref:
    input:
        vcf = rules.align_vcf_to_reference.output.vcf,
        ref = rules.filter_for_eur.output.vcf,
        beagle = rules.download_beagle.output.beagle,
    output:
        vcf = "results/processed_data/lit/imputation/lit.chr{i}.imputed.vcf.gz",
    resources:
        cpus = 4, mem_gb = 32,
    shell:
        """
        java -Xmx4g -jar {input.beagle} gt={input.vcf} ref={input.ref} out=results/processed_data/lit/imputation/lit.chr{wildcards.i}.imputed nthreads={resources.cpus} ne=10000
        tabix -p vcf {output.vcf}
        """


rule merge_lit_yri:
    input:
        vcf1 = rules.beagle_imputation_with_ref.output.vcf,
        vcf2 = rules.filter_for_yri.output.vcf,
    output:
        vcf = "results/processed_data/lit/snps/lit.chr{i}.biallelic.snps.vcf.gz",
    shell:
        """
        bcftools merge {input.vcf1} {input.vcf2} | bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule create_lit_samples:
    input:
        vcf = rules.plink_to_vcf.output.vcf,
        yri_samples = rules.extract_yoruba_eur_samples.output.yri_samples,
        nea_vcf = "resources/data/NeaAltai/not_filtered/chr21_mq25_mapab100.vcf.gz",
        den_vcf = "resources/data/Denisova/not_filtered/chr21_mq25_mapab100.vcf.gz",
    output:
        ref = "results/processed_data/lit/samples/lit.ref.samples.txt",
        tgt = "results/processed_data/lit/samples/lit.tgt.samples.txt",
        nea = "results/processed_data/lit/samples/lit.nea.samples.txt",
        den = "results/processed_data/lit/samples/lit.den.samples.txt",
    shell:
        """
        awk '{{print "YRI\\t"$0}}' {input.yri_samples} > {output.ref}
        bcftools query -l {input.vcf} | awk '{{print "LIT\\t"$0}}' > {output.tgt}
        bcftools query -l {input.nea_vcf} | awk '{{print "NEA\\t"$0}}' > {output.nea}
        bcftools query -l {input.den_vcf} | awk '{{print "DEN\\t"$0}}' > {output.den}
        """


rule add_1src_lit:
    input:
        vcf = "results/processed_data/lit/snps/lit.chr{i}.biallelic.snps.vcf.gz", 
        nea_vcf = rules.filter_nea_genome.output.vcf,
    output:
        vcf = "results/processed_data/lit/1src/lit.nea.chr{i}.vcf.gz",
        index = "results/processed_data/lit/1src/lit.nea.chr{i}.vcf.gz.tbi",
    shell:
        """
        bcftools merge {input.vcf} {input.nea_vcf} | bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule annotate_1src_lit:
    input:
        vcf = rules.add_1src_lit.output.vcf,
        avsnp150 = rules.download_annovar_db.output.hg19_avsnp150,
        dbnsfp42c = rules.download_annovar_db.output.hg19_dbnsfp42c,
    output:
        avinput = "results/annotated_data/lit/lit.nea.chr{i}.avinput",
        txt = "results/annotated_data/lit/lit.nea.chr{i}.hg19_multianno.txt",
    resources:
        cpus=8, mem_gb=32,
    params:
        output_prefix = "results/annotated_data/lit/lit.nea.chr{i}",
    shell:
        """
        bcftools query -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\n" {input.vcf} > {output.avinput}
        resources/tools/annovar/table_annovar.pl {output.avinput} resources/tools/annovar/humandb/ -buildver hg19 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . --thread {resources.cpus}
        """

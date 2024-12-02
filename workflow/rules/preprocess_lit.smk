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


ruleorder: beagle_imputation_with_ref > extract_chromosome


rule plink_to_vcf:
    input:
        ped = rules.download_lithuanian_genomes.output.ped,
        map = rules.download_lithuanian_genomes.output.map,
    output:
        vcf = "results/processed_data/Lithuanians/lit.vcf.gz",
    shell:
        """
        plink --ped {input.ped} --map {input.map} --recode vcf --out results/processed_data/Lithuanians/lit
        bgzip -c results/processed_data/Lithuanians/lit.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_chromosome:
    input:
        vcf = rules.plink_to_vcf.output.vcf,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr{i}.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -r {wildcards.i} | bgzip -c > {output.vcf}
        """


rule align_vcf_to_reference:
    input:
        vcf = rules.extract_chromosome.output.vcf,
        reference = rules.download_1KG_info.output.ref,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr{i}.aligned.vcf.gz",
    script:
        "../scripts/align_reference.py"


rule tabix_index:
    input:
        vcf = rules.align_vcf_to_reference.output.vcf,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr{i}.aligned.vcf.gz.tbi",
    shell:
        """
        tabix -p vcf {input.vcf}
        """


rule extract_yoruba_eur_samples:
    input:
        samples = rules.download_1KG_info.output.samples,
    output:
        yri_samples = "results/processed_data/1KG/yri_samples.txt",
        eur_samples = "results/processed_data/1KG/european_samples.txt",
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
        vcf = "results/processed_data/1KG/chr{i}.EUR.vcf.gz",
    shell:
        """
        bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """


rule filter_for_yri:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
        samples = rules.extract_yoruba_eur_samples.output.yri_samples,
    output:
        vcf = "results/processed_data/1KG/chr{i}.YRI.vcf.gz",
    shell:
        """
        bcftools view -S {input.samples} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """


rule beagle_imputation_with_ref:
    input:
        vcf = rules.align_vcf_to_reference.output.vcf,
        ref = rules.filter_for_eur.output.vcf,
        flag = rules.download_beagle.output.download_flag,
    output:
        vcf = "results/processed_data/Lithuanians/lit.chr{i}.imputed.vcf.gz",
    resources:
        cpus = 2,
    shell:
        """
        java -Xmx4g -jar resources/tools/beagle/beagle.06Aug24.a91.jar gt={input.vcf} ref={input.ref} out=results/processed_data/Lithuanians/lit.chr{wildcards.i}.imputed nthreads={resources.cpus} ne=10000
        tabix -p vcf {output.vcf}
        """


#rule extract_lit_samples:
#    input:
#        vcf = rules.beagle_imputation_with_ref.output.vcf,
#        yri_samples = rules.extract_yoruba_eur_samples.output.yri_samples,
#    output:
#        yri_samples = "results/processed_data/Lithuanians/yri.samples.txt",
#        lit_samples = "results/processed_data/Lithuanians/lit.samples.txt",
#    shell:
#        """
#        awk '{{print "YRI\\t"$0}}' {input.yri_samples} > {output.yri_samples}
#        bcftools query -l {input.vcf} | awk '{{print "LIT\\t"$0}}' > {output.lit_samples}
#        """


#rule merge_lit_yri:
#    input:
#        vcf1 = rules.beagle_imputation_with_ref.output.vcf,
#        vcf2 = rules.filter_for_yri.output.vcf,
#    output:
#        vcf = "results/processed_data/Lithuanians/merged_lit_yri.chr6.vcf.gz",
#    shell:
#        """
#        bcftools isec -n=2 -p dir {input.vcf1} {input.vcf2}
#        bgzip -c dir/0000.vcf > dir/0000.vcf.gz
#        tabix -p vcf dir/0000.vcf.gz
#        bgzip -c dir/0001.vcf > dir/0001.vcf.gz
#        tabix -p vcf dir/0001.vcf.gz
#        bcftools merge dir/0000.vcf.gz dir/0001.vcf.gz | bcftools view -v snps -m 2 -M 2 | bgzip -c > {output.vcf}
#        tabix -p vcf {output.vcf}
#        rm -r dir
#        """


#rule merge_lit_yri_nea:
#    input:
#        vcf1 = rules.merge_lit_yri.output.vcf,
#        vcf2 = rules.download_nea_genome.output.vcf,
#    output:
#        vcf = "results/processed_data/Lithuanians/merged_lit_yri_nea.chr6.vcf.gz",
#    shell:
#        """
#        bcftools isec -n=2 -p dir {input.vcf1} {input.vcf2}
#        bgzip -c dir/0000.vcf > dir/0000.vcf.gz
#        tabix -p vcf dir/0000.vcf.gz
#        bgzip -c dir/0001.vcf > dir/0001.vcf.gz
#        tabix -p vcf dir/0001.vcf.gz
#        bcftools merge dir/0000.vcf.gz dir/0001.vcf.gz | bcftools view -v snps -m 2 -M 2 | bgzip -c > {output.vcf}
#        rm -r dir
#        """

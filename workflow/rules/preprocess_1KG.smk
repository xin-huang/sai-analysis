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
        vcf = "results/processed_data/1KG/1KG.chr{i}.biallelic.snps.vcf.gz",
    shell:
        """
        bcftools view {input.vcf} -m 2 -M 2 -v snps | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule merge_1KG_Nea:
    input:
        vcf1 = rules.extract_1KG_biallelic_snps.output.vcf, 
        vcf2 = rules.download_nea_genome.output.vcf,
    output:
        vcf = "results/processed_data/1KG/merged_1KG_nea.chr{i}.vcf.gz",
    params:
        dir = "dir_{i}",
    shell:
        """
        bcftools isec -n=2 -p {params.dir} {input.vcf1} {input.vcf2}
        bgzip -c {params.dir}/0000.vcf > {params.dir}/0000.vcf.gz
        tabix -p vcf {params.dir}/0000.vcf.gz
        bgzip -c {params.dir}/0001.vcf > {params.dir}/0001.vcf.gz
        tabix -p vcf {params.dir}/0001.vcf.gz
        bcftools merge {params.dir}/0000.vcf.gz {params.dir}/0001.vcf.gz | bcftools view -v snps -m 2 -M 2 | bcftools annotate -x 'FORMAT' | bcftools annotate -x 'INFO' | bgzip -c > {output.vcf}
        rm -r {params.dir}
        """

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


rule analyze_2src_samples:
    input:
        vcf = rules.add_2src_1KG.output.vcf,
        ref = "results/processed_data/1KG/samples/1KG.ref.samples.txt",
        tgt = "results/processed_data/1KG/samples/1KG.tgt.samples.txt",
        src = "results/processed_data/1KG/samples/1KG.nea_den.samples.txt",
        anc_alleles = rules.fa2bed.output.bed,
    output:
        scores = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{stat}.{allele_type}.scores.tsv",
    params:
        chr_name = "{i}",
        win_len = 40000,
        win_step = 40000,
        w = "{w}",
        y = "{y} {z}",
        anc_alleles = lambda wildcards: "" if wildcards.allele_type != "derived" else f"--anc-alleles results/polarized_data/anc_info/hg19.chr{wildcards.i}.anc.alleles.bed",
    resources:
        mem_gb = 16,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --w {params.w} --y {params.y} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --stat {wildcards.stat} {params.anc_alleles}
        """


rule get_2src_samples_outliers:
    input:
        scores = expand("results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{stat}.{allele_type}.scores.tsv", i=list(range(1,23)), allow_missing=True),
    output:
        all_scores = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{stat}.{allele_type}.scores.tsv",
        outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{stat}.{allele_type}.outliers.tsv",
    params:
        outlier_quantile = 0.99,
    shell:
        """
        set +o pipefail
        cat {input.scores} | head -1 > {output.all_scores}
        cat {input.scores} | grep -v "Chrom" >> {output.all_scores}
        sai outlier --score {output.all_scores} --output {output.outliers} --quantile {params.outlier_quantile}
        """

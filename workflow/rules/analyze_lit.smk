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


rule analyze_lit_samples:
    input:
        vcf = "results/processed_data/lit/1src/lit.nea.chr{i}.vcf.gz",
        ref = "results/processed_data/lit/samples/lit.ref.samples.txt",
        tgt = "results/processed_data/lit/samples/lit.tgt.samples.txt",
        src = "results/processed_data/lit/samples/lit.nea.samples.txt",
        anc_alleles = rules.fa2bed.output.bed,
    output:
        scores = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.chr{i}.w_{w}_y_{y}.{stat}.{allele_type}.scores.tsv",
    params:
        chr_name = "{i}",
        win_len = 50000,
        win_step = 10000,
        w = "{w}",
        y = "{y}",
        anc_alleles = lambda wildcards: "" if wildcards.allele_type != "derived" else f"--anc-alleles results/polarized_data/anc_info/hg19.chr{wildcards.i}.anc.alleles.bed",
    resources:
        mem_gb = 16,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --w {params.w} --y {params.y} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --stat {wildcards.stat} {params.anc_alleles}
        """


rule get_lit_outliers:
    input:
        scores = expand("results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.chr{i}.w_{w}_y_{y}.{stat}.{allele_type}.scores.tsv", i=list(range(1,23)), allow_missing=True),
    output:
        all_scores = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{stat}.{allele_type}.scores.tsv",
        outliers = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{stat}.{allele_type}.outliers.tsv",
    params:
        outlier_quantile = 0.999,
    shell:
        """
        set +o pipefail
        cat {input.scores} | head -1 > {output.all_scores}
        cat {input.scores} | grep -v "Chrom" >> {output.all_scores}
        sai outlier --score {output.all_scores} --output {output.outliers} --quantile {params.outlier_quantile}
        """

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


rule analyze_1KG_1src_samples:
    input:
        vcf = rules.merge_1KG_src.output.vcf,
        ref = rules.create_ref_tgt_samples.output.ref,
        tgt = rules.create_ref_tgt_samples.output.tgt,
        src = rules.extract_src_samples.output.src_samples,
    output:
        scores = "results/sai/1src/1KG/{src}/w_{w}_x_{x}_y_{y}/1KG.{src}.chr{i}.w_{w}_x_{x}_y_{y}.scores.tsv",
    params:
        chr_name = "{i}",
        win_len = 40000,
        win_step = 40000,
        w = "{w}",
        x = "{x}",
        y = "{y}",
        q = 0.95,
    resources:
        cpus = 1, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus}
        """


rule analyze_lit_1src_samples:
    input:
        vcf = rules.merge_lit_yri_src.output.vcf,
        ref = rules.extract_lit_samples.output.yri_samples,
        tgt = rules.extract_lit_samples.output.lit_samples,
        src = rules.extract_src_samples.output.src_samples,
    output:
        scores = "results/sai/1src/lit/{src}/w_{w}_x_{x}_y_{y}/lit.{src}.chr{i}.w_{w}_x_{x}_y_{y}.scores.tsv",
    params:
        chr_name = "{i}",
        win_len = 50000,
        win_step = 10000,
        w = "{w}",
        x = "{x}",
        y = "{y}",
        q = 0.95,
    resources:
        cpus = 1, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus}
        """


rule get_1src_samples_outliers:
    input:
        scores = expand("results/sai/1src/{dataset}/{src}/w_{w}_x_{x}_y_{y}/{dataset}.{src}.chr{i}.w_{w}_x_{x}_y_{y}.scores.tsv", i=list(range(1,23)), allow_missing=True),
    output:
        all_scores = "results/sai/1src/{dataset}/{src}/w_{w}_x_{x}_y_{y}/{dataset}.{src}.w_{w}_x_{x}_y_{y}.scores.tsv",
        u_outliers = "results/sai/1src/{dataset}/{src}/w_{w}_x_{x}_y_{y}/{dataset}.{src}.w_{w}_x_{x}_y_{y}_U_outliers.tsv",
        q_outliers = "results/sai/1src/{dataset}/{src}/w_{w}_x_{x}_y_{y}/{dataset}.{src}.w_{w}_x_{x}_y_{y}_Q95_outliers.tsv",
    params:
        outlier_quantile = 0.99,
    shell:
        """
        set +o pipefail
        cat {input.scores} | head -1 > {output.all_scores}
        cat {input.scores} | grep -v "Chrom" >> {output.all_scores}
        sai outlier --score {output.all_scores} --output {output.u_outliers} --stat-type U --quantile {params.outlier_quantile}
        sai outlier --score {output.all_scores} --output {output.q_outliers} --stat-type Q --quantile {params.outlier_quantile}
        """

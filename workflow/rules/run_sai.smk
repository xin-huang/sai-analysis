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


rule analyze_lit_samples:
    input:
        vcf = rules.merge_lit_yri_nea.output.vcf,
        ref = rules.extract_lit_samples.output.yri_samples,
        tgt = rules.extract_lit_samples.output.lit_samples,
        src = rules.extract_nea_samples.output.samples,
    output:
        scores = "results/sai/Lithuanians/lit.scores.txt",
    params:
        chr_name = 6,
        win_len = 50000,
        win_step = 10000,
        w = 0.3,
        x = 0.5,
        y = 1,
    resources:
        cpus = 32, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus}
        """
